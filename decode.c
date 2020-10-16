#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> // for allocation in main

typedef unsigned int uint;
typedef uint8_t   u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int16_t  s16;
typedef int32_t  s32;
typedef int64_t  s64;

#include "farmdec.h"

static Inst UNKNOWN_INST = {
	op: A64_UNKNOWN,
	// all other fields: zero
};

static Inst errinst(char *err) {
	Inst inst = UNKNOWN_INST;
	inst.op = A64_ERROR,
	inst.error = err;
	return inst;
}

Cond get_cond(u8 flags) {
	return (flags >> 4) & 0b1111;
}

static u8 set_cond(u8 flags, Cond cond) {
	cond &= 0xF;
	flags &= 0x0F;
	return (cond << 4) | flags;
}

static u8 invert_cond(u8 flags) {
	Cond cond = get_cond(flags);
	return set_cond(flags, cond ^ 0b001); // invert LSB
}

// Addressing mode, for Loads and Stores.
AddrMode get_addrmode(u8 flags) {
	return (AddrMode)((flags >> 5) & 0b111);
}

static u8 set_addrmode(u8 flags, AddrMode mode) {
	return ((mode&0b111) << 5) | (flags&0b11111);
}

// How much memory to load/store (access size) and whether to sign-
// or zero-extend the value.
ExtendType get_mem_extend(u8 flags) {
	return (ExtendType)((flags >> 2) & 0b111);
}

static u8 set_mem_extend(u8 flags, ExtendType memext) {
	return ((memext&0b111) << 2) | (flags&0b11100011);
}

// The destination register Rd, if present, occupies bits 0..4.
// Register 31 is treated as the Zero/Discard register ZR/WZR.
static Reg regRd(u32 binst) {
	return binst & 0b11111;
}

// Register 31 is treated as the stack pointer SP.
static Reg regRdSP(u32 binst) {
	Reg rd = binst & 0b11111;
	return (rd == 31) ? STACK_POINTER : rd;
}

// The first operand register Rn, if present, occupies bits 5..9.
// Register 31 is treated as the Zero/Discard register ZR/WZR.
static Reg regRn(u32 binst) {
	return (binst >> 5) & 0b11111;
}

// Register 31 is treated as the stack pointer SP.
static Reg regRnSP(u32 binst) {
	Reg rn = (binst >> 5) & 0b11111;
	return (rn == 31) ? STACK_POINTER : rn;
}

// The second operand register Rm, if present, occupies bits 16..20.
// Register 31 is treated as the Zero/Discard register ZR/WZR.
static Reg regRm(u32 binst) {
	return (binst >> 16) & 0b11111;
}

// Register 31 is treated as the stack pointer SP.
static Reg regRmSP(u32 binst) {
	Reg rm = (binst >> 16) & 0b11111;
	return (rm == 31) ? STACK_POINTER : rm;
}

// sext sign-extends the b-bits number in x to 64 bit. The upper (64-b) bits
// must be zero. Seldom needed, but fiddly.
//
// Taken from https://graphics.stanford.edu/~seander/bithacks.html#VariableSignExtend
static s64 sext(u64 x, u8 b) {
	s64 mask = 1U << (b - 1);
	return (((s64)x) ^ mask) - mask;
}

static Inst find_bfm_alias(Op op, bool w32, Reg rd, Reg rn, u8 immr, u8 imms);

static Inst data_proc_imm(u32 binst) {
	Inst inst = UNKNOWN_INST;

	u32 op01 = (binst >> 22) & 0b1111; // op0 and op1 together
	u32 top3 = (binst >> 29) & 0b111;

	enum {
		Unknown,
		PCRelAddr,
		AddSubTags,
		AddSub,
		Logic,
		Move,
		Bitfield,
		Extract
	} kind = Unknown;

	switch (op01) {
	case 0b0000:
	case 0b0001:
	case 0b0010:
	case 0b0011: // 00xx
		kind = PCRelAddr;
		break;
	case 0b0110:
	case 0b0111: // 011x
		kind = AddSubTags;
		break;
	case 0b0100:
	case 0b0101: // 010x
		kind = AddSub;
		break;
	case 0b1000:
	case 0b1001: // 100x
		kind = Logic;
		break;
	case 0b1010:
	case 0b1011: // 101x
		kind = Move;
		break;
	case 0b1100:
	case 0b1101: // 110x
		kind = Bitfield;
		break;
	case 0b1110:
	case 0b1111: // 111x
		kind = Extract;
		break;
	}

	// Bit 31 (sf) controls length of registers (0 → 32 bit, 1 → 64 bit)
	// for most of these data processing operations.
	if ((top3 & 0b100) == 0)
		inst.flags |= W32;

	switch (kind) {
	case Unknown:
		return UNKNOWN_INST;

	case PCRelAddr: {
		if ((top3 & 0b100) == 0)
			inst.op = A64_ADR;
		else
			inst.op = A64_ADRP;

		inst.flags &= ~W32; // no 32-bit variant of these

		// First, just extract the immediate.
		u32 immhi = binst & (0b111111111111111111 << 5);
		u32 immlo = top3 & 0b011;
		u32 uimm = (immhi >> (5-2)) | immlo; // pos(immhi) = 5; 2: len(immlo)

		u64 scale = (inst.op == A64_ADRP) ? 4096 : 1; // ADRP: Page (4K) Address
		s64 simm = sext(scale * uimm, 21);            // PC-relative → signed
		inst.offset = simm;

		inst.rd = regRd(binst);
		break;
	}
	case AddSubTags:
		return errinst("ADDG, SUBG not supported");

	case AddSub: {
		bool is_add = (top3 & 0b010) == 0;
		inst.op = (is_add) ? A64_ADD_IMM : A64_SUB_IMM;
		if (top3 & 0b001)
			inst.flags |= SET_FLAGS;

		u64 unshifted_imm = (binst >> 10) & 0b111111111111;
		bool shift_by_12 = (binst & (1 << 22)) > 0;
		inst.imm = (shift_by_12) ? unshifted_imm << 12 : unshifted_imm;

		inst.rd = regRdSP(binst);
		inst.rn = regRnSP(binst);

		if (inst.op == A64_ADD_IMM && !shift_by_12 && unshifted_imm == 0 &&
			((inst.rd == STACK_POINTER) || (inst.rn == STACK_POINTER))) {
			inst.op = A64_MOV_SP;
		}
		break;
	}
	case Logic: {
		switch (top3 & 0b011) {
		case 0b00: inst.op = A64_AND_IMM; break;
		case 0b01: inst.op = A64_ORR_IMM; break;
		case 0b10: inst.op = A64_EOR_IMM; break;
		case 0b11:
			inst.op = (regRd(binst) == ZERO_REG) ? A64_TST_IMM : A64_AND_IMM;
			inst.flags |= SET_FLAGS;
			break;
		}

		// There's a bit swap: the bits in the instruction are N, immr, imms,
		// from which the immediate imm = N:imms:immr is created.
		u64 immr = binst & (0b111111 << 16); // low 6 bit
		u64 imms = binst & (0b111111 << 10); // high 6 bit
		u64 N = (inst.w32) ? 0 : binst & (1 << 22); // N is MSB of imm for 64-bit variants
		inst.imm = N >> (22-6-6) | imms >> (10-6) | immr >> 16;

		inst.rd = regRdSP(binst);
		inst.rn = regRn(binst);
		break;
	}
	case Move: {
		switch (top3 & 0b011) { // opc
		case 0b00: inst.op = A64_MOVN; break;
		case 0b01: return UNKNOWN_INST;
		case 0b10: inst.op = A64_MOVZ; break;
		case 0b11: inst.op = A64_MOVK; break;
		}

		u8 hw = (binst >> 21) & 0b11;
		u64 imm16 = (binst >> 5) & 0xFFFF;
		inst.mov_wide.imm16 = imm16;
		inst.mov_wide.lsl = 16 * hw;

		inst.rd = regRd(binst);

		// Unshifted or nonzero immediate... there might be a MOV alias.
		if (imm16 != 0 || hw == 0) {
			switch (inst.op) {
			case A64_MOVN:
				if ((inst.flags | W32) && imm16 == 0xFFFF) {
					break; // for 32-bits variant, imm16 may not be all-ones
				}

				inst.op = A64_MOV_IMM;
				inst.imm = ~ (imm16 << inst.mov_wide.lsl);
				break;
			case A64_MOVZ:
				inst.op = A64_MOV_IMM;
				inst.imm = imm16 << inst.mov_wide.lsl;
				break;
			case A64_MOVK:
				break;
			default:
				// Shut up impossible unhandled enum value warning.
				return errinst("data_proc_imm/Move: cannot happen");
			}
		}

		break;
	}
	case Bitfield: {
		Op op = A64_UNKNOWN;

		switch (top3 & 0b011) { // base instructions
		case 0b00: op = A64_SBFM; break;
		case 0b01: op = A64_BFM;  break;
		case 0b10: op = A64_UBFM; break;
		default:
			return errinst("data_proc_imm/Bitfield: neither SBFM, BFM or UBFM");
		}

		bool w32 = (inst.flags & W32);
		u8 immr = (binst >> 16) & 0b111111;
		u8 imms = (binst >> 10) & 0b111111;
		Reg rd = regRd(binst);
		Reg rn = regRn(binst);
		inst = find_bfm_alias(op, w32, rd, rn, immr, imms);
		break;
	}
	case Extract: {
		inst.op = A64_EXTR;
		inst.imm = (binst >> 10) & 0b111111;
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		if (inst.rn == inst.rm) {
			inst.op = A64_ROR_IMM;
			inst.rm = 0; // unused for ROR_IMM → clear again
		}
		break;
	}
	}

	return inst;
}

// There is no BFM, SBFM or UBFM instruction without a more specific alias.
// Returns a complete instruction, with flags and fields properly set.
static Inst find_bfm_alias(Op op, bool w32, Reg rd, Reg rn, u8 immr, u8 imms) {
	Inst inst = UNKNOWN_INST;
	u8 all_ones = (w32) ? 31 : 63; // for ASR, LSL, LSR
	u8 bits = (w32) ? 32 : 64;     // for BFC, BFI, LSL, SBFIZ, UBFIZ

	inst.rd = rd;
	inst.rn = rn;
	inst.flags |= W32;

	// General note:
	//
	//   immr = -lsb MOD bits
	// ≡ immr = bits - lsb   (∀lsb <= bits)
	// ≡ lsb  = bits - immr

	if (op == A64_BFM) {
		// BFM Rd, Rn, immr=#lsb, imms=#(lsb+width-1) → BFXIL Rd, Rn, #lsb, #width
		if (imms >= immr) {
			inst.op = A64_BFXIL;
			inst.bfm.lsb = immr;
			inst.bfm.width = imms - immr + 1;
			return inst;
		}

		// BFM Rd, ZR, immr=#(-lsb MOD bits), imms=#(width-1) → BFC Rd, #lsb, #width
		// BFM Rd, Rn, immr=#(-lsb MOD bits), imms=#(width-1) → BFI Rd, Rn, #lsb, #width
		inst.op = (rn == ZERO_REG) ? A64_BFC : A64_BFI;
		inst.bfm.lsb = bits - immr;
		inst.bfm.width = imms + 1;
		return inst;
	}

	// SBFM and UBFM are extremely similar except for their sign/opcode.
	bool sign = (op == A64_SBFM);

	// UBFM Rd, Rn, immr=#(-<shift> MOD 32/64), imms=#(31/63 - <shift>) → LSL Rd, Rn, #shift
	if (!sign && imms + 1 == immr && imms != all_ones) {
		inst.op = A64_LSL_IMM;
		inst.imm = all_ones - imms; // 31 - imms or 63 - imms
		return inst;
	}

	// UBFM Rd, Rn, immr=#shift, imms=#31/#63 → LSR Rd, Rn, #shift
	// SBFM Rd, Rn, immr=#shift, imms=#31/#63 → ASR Rd, Rn, #shift
	if (imms == all_ones) {
		inst.op = (sign) ? A64_ASR_IMM : A64_LSR_IMM;
		inst.imm = immr;
		return inst;
	}

	// UBFM Rd, ZR, immr=#(-lsb MOD bits), imms=#(width-1) → UBFIZ Rd, #lsb, #width
	// SBFM Rd, ZR, immr=#(-lsb MOD bits), imms=#(width-1) → SBFIZ Rd, #lsb, #width
	if (imms < immr) {
		inst.op = (sign) ? A64_SBFIZ : A64_UBFIZ;
		inst.bfm.lsb = bits - immr;
		inst.bfm.width = imms + 1;
		return inst;
	}

	// UBFM Rd, Rn, immr=#0, imms=#7  → UXTB Rd, Rn
	// UBFM Rd, Rn, immr=#0, imms=#15 → UXTH Rd, Rn
	// SBFM Rd, Rn, immr=#0, imms=#7  → SXTB Rd, Rn
	// SBFM Rd, Rn, immr=#0, imms=#15 → SXTH Rd, Rn
	// SBFM Rd, Rn, immr=#0, imms=#31 → SXTW Rd, Rn
	if (immr == 0) {
		switch (imms) {
		case  7: inst.op = A64_EXTEND; inst.extend.type = (sign) ? SXTB : UXTB; return inst;
		case 15: inst.op = A64_EXTEND; inst.extend.type = (sign) ? SXTH : UXTH; return inst;
		case 31:
			inst.op = A64_EXTEND;
			inst.extend.type = SXTW;
			return (sign) ? inst : errinst("find_bfm_alias: there is no UXTW instruction");
		}
	}

	// This catches all the cases for which no more specific alias exists.
	//
	// UBFM Rd, Rn, immr=#lsb, imms=#(lsb+width-1) → UBFX Rd, Rn, #lsb, #width
	// SBFM Rd, Rn, immr=#lsb, imms=#(lsb+width-1) → SBFX Rd, Rn, #lsb, #width
	inst.op = (sign) ? A64_SBFX : A64_UBFX;
	inst.bfm.lsb = immr;
	inst.bfm.width = imms - immr + 1;
	return inst;
}

static Inst branches(u32 binst) {
	Inst inst = UNKNOWN_INST;
	u32 top7 = (binst >> 25) & 0b1111111;

	enum {
		Unknown,
		CondBranch,      // B.cond
		System,          // many instructions, mostly irrelevant
		UncondBranchReg, // BR, BRL, RET
		UncondBranch,    // B, BL
		CmpAndBranch,    // CBZ, CBNZ
		TestAndBranch    // TBZ, TBNZ
	} kind = Unknown;

	// XXX there are more patterns for less important instructions.
	// XXX order everything by the table in the documentation.
	switch (top7) {
	case 0b0101010:
		kind = CondBranch;
		break;

	case 0b1101010:
		// Exception generation, hints, barriers, PSTATE, sys instrs, sys reg moves.
		// Most of them are entirely unimportant for us, the probable exceptions being
		// barriers and CFINV.
		kind = System;
		break;

	case 0b1101011:
		kind = UncondBranchReg;
		break;

	case 0b0001010:
	case 0b0001011:
	case 0b1001010:
	case 0b1001011: // x00101x
		kind = UncondBranch;
		break;

	case 0b0011010:
	case 0b1011010: // x011010
		kind = CmpAndBranch;
		break;

	case 0b0011011:
	case 0b1011011: // x011011
		kind = TestAndBranch;
		break;

	default:
		return UNKNOWN_INST;
	}

	// Regarding "immX * 4": branch immediates are word offsets, but we want
	// byte offsets, so shift twice to the left.

	switch (kind) {
	case Unknown:
		return UNKNOWN_INST;

	case CondBranch:
		inst.op = A64_BCOND;
		inst.flags = set_cond(inst.flags, binst & 0b1111);
		inst.offset = sext(((binst >> 5) & 0b1111111111111111111) << 2, 19+2); // imm19 * 4
		break;

	case System: {
		bool exc = !((binst >> 24) & 1); // exception → zero bit

		if (exc) {
			u8 opcLL = (((binst >> 21) & 0b111) << 2) | (binst & 0b11); // opc(3):LL(2)
			switch (opcLL) {
			case 0b00001: inst.op = A64_SVC; break; // system call
			case 0b00010: inst.op = A64_HVC; break;
			case 0b00011: inst.op = A64_SMC; break;
			case 0b00100: inst.op = A64_BRK; break;
			case 0b01000: inst.op = A64_HLT; break;
			case 0b10101: inst.op = A64_DCPS1; break;
			case 0b10110: inst.op = A64_DCPS2; break;
			case 0b10111: inst.op = A64_DCPS3; break;
			default:
				return errinst("branches/System/exceptions: bad instruction");
			}

			inst.imm = (binst >> 5) & 0xFFFF; // imm16
			break;
		}
		return errinst("System instructions not supported"); // XXX _some_ should be decoded. barriers, CFINV
	}
	case UncondBranchReg: {
		u32 op = (binst >> 21) & 0b11;
		switch (op) {
		case 0b00: inst.op = A64_BR;  break;
		case 0b01: inst.op = A64_BLR; break;
		case 0b10: inst.op = A64_RET; break;
		default:
			return UNKNOWN_INST;
		}

		inst.rn = regRn(binst);
		break;
	}
	case UncondBranch:
		inst.op = (binst & (1<<31)) ? A64_BL : A64_B; // MSB = 1 → BL
		inst.offset = sext((binst & 0b11111111111111111111111111) << 2, 26+2); // imm26 * 4
		break;

	case CmpAndBranch: {
		if ((top7 & 0b1000000) == 0)
			inst.flags |= W32;

		bool zero = (binst & (1 << 24)) == 0;
		inst.op = (zero) ? A64_CBZ : A64_CBNZ;

		inst.offset = 4 * sext(((binst >> 5) & 0b1111111111111111111) << 2, 19+2); // imm19 * 4
		inst.rn = binst & 0b11111; // Rt; not modified → Inst.rn, not .rd
		break;
	}
	case TestAndBranch: {
		if ((top7 & 0b1000000) == 0)
			inst.flags |= W32;

		bool zero = (binst & (1 << 24)) == 0;
		inst.op = (zero) ? A64_TBZ : A64_TBNZ;

		inst.tbz.offset = sext(((binst >> 5) & 0b11111111111111) << 2, 14+2); // imm14 * 4
		u32 b40 = (binst >> 19) & 0b11111;
		u32 b5 = binst & (1<<31);
		inst.tbz.bit = (b5 >> (31-5)) | b40; // b5:b40
		inst.rn = binst & 0b11111;           // Rt; not modified → Inst.rn, not .rd
		break;
	}
	}

	return inst;
}

static Inst data_proc_reg(u32 binst) {
	Inst inst = UNKNOWN_INST;
	u32 top3 = (binst >> 29) & 0b111;
	u32 op = (binst >> 21) & 0b11111111; // op1:101:op2 bits 21..27

	enum {
		Unknown,
		DataProc2,
		DataProc1,
		LogicalShifted,
		AddSubShifted,
		AddSubExtended,
		AddSubCarryAndFlagManipulation,
		CondCompare,
		CondSelect,
		DataProc3
	} kind = Unknown;

	switch (op) {
	case 0b11010110: // dp2, dp1
		kind = (top3 & 0b010) ? DataProc1 : DataProc2;
		break;
	case 0b01010000:
	case 0b01010001:
	case 0b01010010:
	case 0b01010011:
	case 0b01010100:
	case 0b01010101:
	case 0b01010110:
	case 0b01010111: // 01010xxx
		kind = LogicalShifted;
		break;
	case 0b01011000:
	case 0b01011010:
	case 0b01011100:
	case 0b01011110: // 01011xx0
		kind = AddSubShifted;
		break;
	case 0b01011001:
	case 0b01011011:
	case 0b01011101:
	case 0b01011111: // 01011xx0
		kind = AddSubExtended;
		break;
	case 0b11010000:
		// ADC, ADCS, SBC, SBCS, RMIF, SETF8, SETF16
		// Really similar bit patterns, so treat it as one kind.
		kind = AddSubCarryAndFlagManipulation;
		break;
	case 0b11010010:
		kind = CondCompare;
		break;
	case 0b11010100:
		kind = CondSelect;
		break;
	case 0b11011000:
	case 0b11011001:
	case 0b11011010:
	case 0b11011011:
	case 0b11011100:
	case 0b11011101:
	case 0b11011110:
	case 0b11011111: // 11011xxx
		kind = DataProc3;
		break;
	default:
		return UNKNOWN_INST;
	}

	// Bit 31 (sf) controls length of registers (0 → 32 bit, 1 → 64 bit)
	// for _all_ register data processing instructions. Instructions that
	// have no 32-bit variant still have the bit, but it's always set.
	if ((top3 & 0b100) == 0)
		inst.flags |= W32;

	switch (kind) {
	case Unknown:
		return UNKNOWN_INST;

	case DataProc2: {
		switch ((binst >> 10) & 0b111111) {
		case 0b000000:
			inst.op = A64_SUBP;
			if (top3 & 0b001)
				inst.flags |= SET_FLAGS; // SUBPS
			inst.rd = regRd(binst);
			inst.rn = regRnSP(binst);
			inst.rm = regRmSP(binst);
			return inst;
		case 0b000010: inst.op = A64_UDIV; break;
		case 0b000011: inst.op = A64_SDIV; break;
		case 0b001000: inst.op = A64_LSLV; break; // lowest two bits → enum Shift
		case 0b001001: inst.op = A64_LSRV; break; // --
		case 0b001010: inst.op = A64_ASRV; break; // --
		case 0b001011: inst.op = A64_RORV; break; // --
		case 0b010000: inst.op = A64_CRC32B; break; // lowest two bits → enum Size
		case 0b010001: inst.op = A64_CRC32H; break;
		case 0b010010: inst.op = A64_CRC32W; break;
		case 0b010011: inst.op = A64_CRC32X; break;
		case 0b010100: inst.op = A64_CRC32CB; break;
		case 0b010101: inst.op = A64_CRC32CH; break;
		case 0b010110: inst.op = A64_CRC32CW; break;
		case 0b010111: inst.op = A64_CRC32CX; break;
		}

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		break;
	}
	case DataProc1:
		switch ((binst >> 10) & 0b1111111111) {
		case 0b000: inst.op = A64_RBIT; break;
		case 0b001: inst.op = A64_REV16; break;
		case 0b010: inst.op = (inst.flags & W32) ? A64_REV : A64_REV32; break;
		case 0b100: inst.op = A64_CLZ; break;
		case 0b101: inst.op = A64_CLS; break;
		}

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;

	case LogicalShifted: {
		bool negate = (binst >> 21) & 1;
		switch (top3 & 0b011) {
		case 0b00: inst.op = (negate) ? A64_BIC : A64_AND_SHIFTED; break;
		case 0b01: inst.op = (negate) ? A64_ORN : A64_ORR_SHIFTED; break;
		case 0b10: inst.op = (negate) ? A64_EON : A64_EOR_SHIFTED; break;
		case 0b11:
			inst.op = (negate) ? A64_BIC : A64_AND_SHIFTED;
			inst.flags |= SET_FLAGS; // ANDS and BICS represented using SET_FLAGS bit
			break;
		}

		u32 shift = (binst >> 22) & 0b11;
		u32 imm6 = (binst >> 10) & 0b111111;
		inst.shift.type = shift;
		inst.shift.amount = imm6;

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);

		if (inst.op == A64_AND_SHIFTED && (inst.flags & SET_FLAGS) && inst.rd == ZERO_REG)
			inst.op = A64_TST_SHIFTED;
		else if (inst.op == A64_ORR_SHIFTED && shift == 0 && imm6 == 0 && inst.rn == ZERO_REG)
			inst.op = A64_MOV_REG;
		else if (inst.op == A64_ORN && inst.rn == ZERO_REG)
			inst.op = A64_MVN;
		break;
	}
	case AddSubShifted:
		if (top3 & 0b001)
			inst.flags |= SET_FLAGS;

		switch (top3 & 0b011) {
		case 0b00:
		case 0b01:
			inst.op = A64_ADD_SHIFTED;
			break;
		case 0b10:
		case 0b11:
			inst.op = A64_SUB_SHIFTED;
			break;
		}

		inst.shift.type = (binst >> 22) & 0b11;
		inst.shift.amount = (binst >> 10) & 0b111111;

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);

		if (inst.rd == ZERO_REG && (inst.flags & SET_FLAGS)) {
			switch (inst.op) {
			case A64_ADD_SHIFTED: inst.op = A64_CMN_SHIFTED; break;
			case A64_SUB_SHIFTED: inst.op = A64_CMP_SHIFTED; break;
			default:
				break; // impossible; shuts up warning
			}
		} else if (inst.op == A64_SUB_SHIFTED && inst.rn == ZERO_REG) {
			inst.op = A64_NEG;
		}
		break;

	case AddSubExtended: {
		if (top3 & 0b001)
			inst.flags |= SET_FLAGS;

		switch (top3 & 0b011) {
		case 0b00:
		case 0b01:
			inst.op = A64_ADD_EXT;
			break;
		case 0b10:
		case 0b11:
			inst.op = A64_SUB_EXT;
			break;
		}

		inst.extend.type = (binst >> 13) & 0b111; // three bits: sign(1):size(2)
		inst.extend.lsl = (binst >> 10) & 0b111;  // optional LSL amount

		// Unlike AddSubShifted, R31 is intepreted as the Stack Pointer.
		inst.rd = regRdSP(binst);
		inst.rn = regRnSP(binst);
		inst.rm = regRm(binst);

		if (inst.rd == ZERO_REG && (inst.flags & SET_FLAGS)) {
			switch (inst.op) {
			case A64_ADD_EXT: inst.op = A64_CMN_EXT; break;
			case A64_SUB_EXT: inst.op = A64_CMP_EXT; break;
			default:
				break; // impossible; shuts up warning
			}
		}
		break;
	}
	case AddSubCarryAndFlagManipulation: {
		if (top3 & 0b001)
			inst.flags |= SET_FLAGS;

		switch ((binst >> 9) & 0b111111) {
		case 0b000000: {
			bool sub = top3 & 0b010;
			inst.op = (sub) ? A64_SBC : A64_ADC;
			inst.rd = regRd(binst);
			inst.rn = regRn(binst);
			inst.rm = regRm(binst);
			if (inst.op == A64_SBC && inst.rn == ZERO_REG) {
				inst.op = A64_NGC;
			}
			break;
		}
		case 0b000001:
		case 0b100001: // x00001
			inst.op = A64_RMIF;
			inst.rn = regRn(binst);
			inst.rmif.mask = binst & 0b1111;
			inst.rmif.ror = (binst >> 15) & 0b111111;
			break;
		case 0b000010:
			inst.op = A64_SETF8;
			inst.rn = regRn(binst);
			break;
		case 0b010010:
			inst.op = A64_SETF16;
			inst.rn = regRn(binst);
			break;
		default:
			return UNKNOWN_INST;
		}

		break;
	}
	case CondCompare: {
		bool positive = (top3 & 0b010) == 0;
		bool immediate = (binst >> 11) & 1;
		if (positive) {
			inst.op = (immediate) ? A64_CCMP_IMM : A64_CCMP_REG;
		} else {
			inst.op = (immediate) ? A64_CCMN_IMM : A64_CCMN_REG;
		}

		inst.flags |= SET_FLAGS;
		inst.flags = set_cond(inst.flags, (binst >> 12) & 0b1111);
		inst.ccmp.nzcv = binst & 0b1111;
		inst.ccmp.imm5 = (immediate) ? (binst >> 16) & 0b11111 : 0;
		break;
	}
	case CondSelect: {
		// Combine the op bit (30) and the op2 field (11-10) to fully
		// determine the selection opcode.
		u32 op = (((binst >> 30) & 1) << 1) | ((binst >> 10) & 0b11);
		switch (op) {
		case 0b000: inst.op = A64_CSEL; break;
		case 0b001: inst.op = A64_CSINC; break;
		case 0b100: inst.op = A64_CSINV; break;
		case 0b101: inst.op = A64_CSNEG; break;
		default:
			return UNKNOWN_INST;
		}

		inst.flags = set_cond(inst.flags, (binst >> 12) & 0b1111);
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);

		if (inst.rm == ZERO_REG && inst.rn == ZERO_REG) {
			switch (inst.op) {
			case A64_CSINC: inst.op = A64_CSET; inst.flags = invert_cond(inst.flags); break;
			case A64_CSINV: inst.op = A64_CSETM; inst.flags = invert_cond(inst.flags); break;
			default:
				break;
			}
		} else if (inst.rm == inst.rn && inst.op == A64_CSNEG) {
			inst.op = A64_CNEG;
			inst.flags = invert_cond(inst.flags);
		} else if (inst.rm == inst.rn && inst.op == A64_CSINC) {
			inst.op = A64_CINC;
			inst.flags = invert_cond(inst.flags);
		}
		break;
	}
	case DataProc3: {
		bool sub = (binst >> 15) & 1;
		switch ((binst >> 21) & 0b111) {
		case 0b000: inst.op = (sub) ? A64_MADD : A64_MSUB; break;
		case 0b001: inst.op = (sub) ? A64_SMADDL : A64_SMSUBL; break;
		case 0b010: inst.op = A64_SMULH; break;
		case 0b101: inst.op = (sub) ? A64_UMADDL : A64_UMSUBL; break;
		case 0b110: inst.op = A64_UMULH; break;
		default:
			return errinst("data_proc_reg/DataProc3: invalid opcode");
		}

		// Only MADD und MSUB have 32-bit variants.
		if ((top3 & 0b100) == 0 && inst.op != A64_MADD && inst.op != A64_MSUB)
			return errinst("data_proc_reg/DataProc3: no 32-bit variant except for MADD and MSUB");

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		inst.ra = (binst >> 10) & 0b11111;

		if (inst.ra == ZERO_REG) {
			switch (inst.op) {
			case A64_MADD:   inst.op =    A64_MUL; break;
			case A64_MSUB:   inst.op =   A64_MNEG; break;
			case A64_SMADDL: inst.op =  A64_SMULL; break;
			case A64_SMSUBL: inst.op = A64_SMNEGL; break;
			case A64_UMADDL: inst.op =  A64_UMULL; break;
			case A64_UMSUBL: inst.op = A64_UMNEGL; break;

			case A64_SMULH:
			case A64_UMULH:
			default:
				break;
			}
		}

		break;
	}
	}

	return inst;
}

static Inst loads_and_stores(u32 binst) {
	Inst inst = UNKNOWN_INST;

	u8 top2 = (binst >> 30) & 0b11;

	// For loads and stores, there is a lot of shared structure,
	// so we don't quite follow the structure of the document like
	// in the other instruction groups.
	enum {
		Unknown,
		SIMDMult,       // LDx, STx
		SIMDSingle,     // --------, LDxR
		Tags,           // LDG, STG -- omitted
		Exclusive,      // Exclusive, Load-acquire/store-release, swaps
		Literal,        // Load from label
		PairNoAlloc,    // LDNP, STNP -- normal access with !nontemporal hint
		Pair,           // LDP, STP
		UnscaledImm,    // LDUR, STUR -- Unscaled signed Immediate
		Register,       // LDR, STR
		Unprivileged,   // LDTR, STTR -- omitted
		Atomic,
		PAC             // Pointer Auth -- omitted
	} kind = Unknown;

	// We find the kind as usual, but we also set the addrmode here,
	// so that the switch below can handle the common code regardless
	// of addrmode.
	switch ((binst >> 23) & 0b1111111) {
	case 0b0011000: // LDx, STx mult. struct.
		kind = SIMDMult;
		inst.flags = set_addrmode(inst.flags, AM_SIMPLE);
		break;
	case 0b0011001: // ----- (post-indexed)
		kind = SIMDMult;
		inst.flags = set_addrmode(inst.flags, AM_POST);
		break;
	case 0b0011010: // LDx, STx single struct.
		kind = SIMDSingle;
		inst.flags = set_addrmode(inst.flags, AM_SIMPLE);
		break;
	case 0b0011011: // ----- (post-indexed)
		kind = SIMDSingle;
		inst.flags = set_addrmode(inst.flags, AM_POST);
		break;
	case 0b0110010:
	case 0b0110011: // 011001x -- memory tags (UNIMPLEMENTED)
		kind = Tags;
		break;
	case 0b0010000:
	case 0b0010001: // 001000x -- exclusive
		kind = Exclusive;
		break;
	case 0b0110000:
	case 0b0110001:
	case 0b0111000:
	case 0b0111001: // 011x00x -- load literal
		kind = Literal;
		inst.flags = set_addrmode(inst.flags, AM_LITERAL);
		break;
	case 0b1010000:
	case 0b1011000: // 101x000 -- no-allocate pair
		kind = PairNoAlloc;
		inst.flags = set_addrmode(inst.flags, AM_OFF_IMM);
		break;
	case 0b1010001:
	case 0b1011001: // 101x001 -- pair (post-indexed)
		kind = Pair;
		inst.flags = set_addrmode(inst.flags, AM_POST);
		break;
	case 0b1010010:
	case 0b1011010: // 101x010 -- pair (signed offset)
		kind = Pair;
		inst.flags = set_addrmode(inst.flags, AM_OFF_IMM);
		break;
	case 0b1010011:
	case 0b1011011: // 101x011 -- pair (pre-indexed)
		kind = Pair;
		inst.flags = set_addrmode(inst.flags, AM_PRE);
		break;
	case 0b1110000:
	case 0b1110001:
	case 0b1111000:
	case 0b1111001: { // 111x00x -- many things, mostly Register
		u32 op4 = (binst >> 10) & 0b11;
		u32 b21 = (binst >> 21) & 1;

		if (b21) {
			switch (op4) {
			case 0b00:
				kind = Atomic;
				inst.flags = set_addrmode(inst.flags, AM_SIMPLE);
				break;
			case 0b10: {
				kind = Register;

				// UXTX (0b011) → no extension → shifted register offset
				u8 option = (binst >> 13) & 0b111;
				inst.flags = set_addrmode(inst.flags, (option == UXTX) ? AM_OFF_REG : AM_OFF_EXT);
				break;
			}
			case 0b01:
			case 0b11:
				kind = PAC;
				break;
			}
			break;
		}

		kind = Register;
		switch (op4) {
		case 0b00:
			kind = UnscaledImm; // LDUR, STUR
			inst.flags = set_addrmode(inst.flags, AM_OFF_IMM);
			break;
		case 0b01: inst.flags = set_addrmode(inst.flags, AM_POST); break;
		case 0b10: kind = Unprivileged;                            break;
		case 0b11: inst.flags = set_addrmode(inst.flags, AM_PRE);  break;
		}
		break;
	}
	case 0b1110010:
	case 0b1110011:
	case 0b1111010:
	case 0b1111011: { // 111x01x -- unsigned immediate
		kind = Register;
		// AM_OFF_IMM: unambiguous since UnscaledImm is a different kind and has other opcodes.
		inst.flags = set_addrmode(inst.flags, AM_OFF_IMM);
		break;
	}
	default:
		return errinst("loads_and_stores: invalid instruction");
	}

	switch (kind) {
	case Unknown:
		return UNKNOWN_INST;
	case SIMDMult: {
		u8 Q = (binst >> 30) & 1;
		bool load = (binst >> 22) & 1;
		VectorArrangement va = (((binst >> 10) & 0b11) << 1) | Q; // size(2):Q(1)
		uint nreg = 0;

		switch ((binst >> 12) & 0b1111) { // opcode(4)
		case 0b0000: inst.op = (load) ? A64_LD4_MULT : A64_ST4_MULT; nreg = 4; break;
		case 0b0010: inst.op = (load) ? A64_LD1_MULT : A64_ST1_MULT; nreg = 4; break;
		case 0b0100: inst.op = (load) ? A64_LD3_MULT : A64_ST3_MULT; nreg = 3; break;
		case 0b0110: inst.op = (load) ? A64_LD1_MULT : A64_ST1_MULT; nreg = 3; break;
		case 0b0111: inst.op = (load) ? A64_LD1_MULT : A64_ST1_MULT; nreg = 1; break;
		case 0b1000: inst.op = (load) ? A64_LD2_MULT : A64_ST2_MULT; nreg = 2; break;
		case 0b1010: inst.op = (load) ? A64_LD1_MULT : A64_ST1_MULT; nreg = 2; break;
		}

		inst.flags = set_mem_extend(inst.flags, va);
		inst.simd_ldst.nreg = nreg;
		inst.ldst.rt = regRd(binst);
		inst.ldst.rn = regRnSP(binst);

		// Post-indexed addrmode, either register-based or an immediate that
		// depends on the number of bytes read/written.
		if (get_addrmode(inst.flags) == AM_POST) {
			inst.ldst.rm = regRm(binst);
			if (inst.ldst.rm == ZERO_REG) {
				switch (Q) {
				case 0: inst.simd_ldst.offset = nreg *  8; break; //  64-bit half vectors
				case 1: inst.simd_ldst.offset = nreg * 16; break; // 128-bit full vectors
				}
			}
		}
		break;
	}
	case SIMDSingle: {
		u8 Q = (binst >> 30) & 1;
		bool load = (binst >> 22) & 1;
		bool even = (binst >> 21) & 1; // R(1): determines evenness of inst: R=0 → LD1, LD3; R=1 → LD2, LD4
		FPSize size = 0;
		u8 size_field = (binst >> 10) & 0b11;         // size(2) is called "size", but is actually part of index
		VectorArrangement va = (size_field << 1) | Q; // size(2):Q(1)
		uint nreg = 0;

		if (load) {
			switch ((binst >> 13) & 0b111) { // opcode(3)
			case 0b000: inst.op = (even) ? A64_LD2_SINGLE : A64_LD1_SINGLE; size = FSZ_B; break;
			case 0b001: inst.op = (even) ? A64_LD4_SINGLE : A64_LD3_SINGLE; size = FSZ_B; break;
			case 0b010: inst.op = (even) ? A64_LD2_SINGLE : A64_LD1_SINGLE; size = FSZ_H; break;
			case 0b011: inst.op = (even) ? A64_LD4_SINGLE : A64_LD3_SINGLE; size = FSZ_H; break;
			case 0b100:
				inst.op = (even) ? A64_LD2_SINGLE : A64_LD1_SINGLE;
				size = (size_field == 0b01) ? FSZ_D : FSZ_S;
				break;
			case 0b101:
				inst.op = (even) ? A64_LD4_SINGLE : A64_LD3_SINGLE;
				size = (size_field == 0b01) ? FSZ_D : FSZ_S;
				break;
			case 0b110: inst.op = (even) ? A64_LD2R       : A64_LD1R; size = (binst >> 10) & 0b11; break;
			case 0b111: inst.op = (even) ? A64_LD4R       : A64_LD3R; size = (binst >> 10) & 0b11; break;
			}
		} else {
			switch ((binst >> 13) & 0b111) { // opcode(3)
			case 0b000: inst.op = (even) ? A64_ST2_SINGLE : A64_ST1_SINGLE; size = FSZ_B; break;
			case 0b001: inst.op = (even) ? A64_ST4_SINGLE : A64_ST3_SINGLE; size = FSZ_B; break;
			case 0b010: inst.op = (even) ? A64_ST2_SINGLE : A64_ST1_SINGLE; size = FSZ_H; break;
			case 0b011: inst.op = (even) ? A64_ST4_SINGLE : A64_ST3_SINGLE; size = FSZ_H; break;
			case 0b100:
				inst.op = (even) ? A64_ST2_SINGLE : A64_ST1_SINGLE;
				size = (size_field == 0b01) ? FSZ_D : FSZ_S;
				break;
			case 0b101:
				inst.op = (even) ? A64_ST4_SINGLE : A64_ST3_SINGLE;
				size = (size_field == 0b01) ? FSZ_D : FSZ_S;
				break;
			// No STxR instructions, since "Store and Replicate to all Lanes" is nonsensical.
			}
		}

		switch (inst.op) {
		case A64_LD1R: case A64_LD1_SINGLE: case A64_ST1_SINGLE: nreg = 1; break;
		case A64_LD2R: case A64_LD2_SINGLE: case A64_ST2_SINGLE: nreg = 2; break;
		case A64_LD3R: case A64_LD3_SINGLE: case A64_ST3_SINGLE: nreg = 3; break;
		case A64_LD4R: case A64_LD4_SINGLE: case A64_ST4_SINGLE: nreg = 4; break;
		default: break;
		}

		// The index for the vector element is tightly encoded using various separate bits.
		// Larger elements imply fewer distinct indexes.
		u8 S = (binst >> 12) & 1;
		switch (size) {
		case FSZ_B: inst.simd_ldst.index = (Q << 3) | (S << 2) | size_field; break;        // Q:S:size    0..15
		case FSZ_H: inst.simd_ldst.index = (Q << 2) | (S << 1) | (size_field >> 1); break; // Q:S:size<1> 0..7
		case FSZ_S: inst.simd_ldst.index = (Q << 1) | S; break;                            // Q:S         0..3
		case FSZ_D: inst.simd_ldst.index = Q; break;                                       // Q           0..1
		case FSZ_Q: break; // impossible
		}

		inst.flags = set_mem_extend(inst.flags, va);
		inst.simd_ldst.nreg = nreg;
		inst.ldst.rt = regRd(binst);
		inst.ldst.rn = regRnSP(binst);

		// Post-indexed addrmode, either register-based or an immediate that
		// depends on the number of bytes read/written.
		if (get_addrmode(inst.flags) == AM_POST) {
			inst.ldst.rm = regRm(binst);
			if (inst.ldst.rm == ZERO_REG) {
				switch (size) {
				case FSZ_B: inst.simd_ldst.offset = nreg * 1; break;
				case FSZ_H: inst.simd_ldst.offset = nreg * 2; break;
				case FSZ_S: inst.simd_ldst.offset = nreg * 4; break;
				case FSZ_D: inst.simd_ldst.offset = nreg * 8; break;
				case FSZ_Q: break; // impossible
				}
			}
		}
		break;
	}
	case Tags:
		return errinst("loads_and_stores: Tag instructions not supported");
	case Exclusive: {
		u8 size = top2;
		bool o0 = (binst >> 15) & 1;
		bool o1 = (binst >> 21) & 1;
		bool load = (binst >> 22) & 1;
		bool o2 = (binst >> 23) & 1;

		if (o1) { // CAS, CASL, CASA, CASAL, CASP, CASPL, CASPA, CASPAL, LDXP, LDAXP, STXP, STLXP
			if (o2) { // CAS*
				inst.op = A64_CAS;
				inst.ldst_order.load = (load) ? MO_ACQUIRE : MO_NONE;
				inst.ldst_order.store = (o0) ? MO_RELEASE : MO_NONE;
				goto common;
			}

			// The size field is used to differentiate CASP and pair instructions.
			// Neither CASP nor the pair instructions allow byte or halfword sizes.
			if (size == SZ_B || size == SZ_H) { // CASP
				inst.op = A64_CASP;
				inst.ldst_order.load = (load) ? MO_ACQUIRE : MO_NONE;
				inst.ldst_order.store = (o0) ? MO_RELEASE : MO_NONE;

				// size+2: size=00 → SZ_W (10), size=01 → SZ_X (11) here.
				size += 2;
				goto common;
			}

			// Exclusive Load/Store Pair
			if (load) {
				inst.op = A64_LDXP;
				inst.ldst_order.load = (o0) ? MO_ACQUIRE : MO_NONE;
			} else {
				inst.op = A64_STXP;
				inst.ldst_order.store = (o0) ? MO_RELEASE : MO_NONE;
				inst.ldst_order.rs = regRm(binst); // status register
			}
			inst.ldst.rt2 = (binst >> 10) & 0b11111;
			goto common;
		}

		// LDXR, LDAXR, STXR, STLXR, STLR, STLLR, LDAR, LDLAR
		bool exclusive = !o2; // LDXR, STXR, ...?

		// The ordered, but non-exclusive LDAR, STLR, ... are represented
		// as LDR/STR with addrmode=AM_SIMPLE and the proper ldst_order
		// values.

		if (load) {
			inst.op = (exclusive) ? A64_LDXR : A64_LDR;
			if (exclusive)
				inst.ldst_order.load = (o0) ? MO_ACQUIRE : MO_NONE;
			else
				inst.ldst_order.load = (o0) ? MO_LO_ACQUIRE : MO_ACQUIRE;
		} else {
			inst.op = (exclusive) ? A64_STXR : A64_STR;
			if (exclusive)
				inst.ldst_order.store = (o0) ? MO_RELEASE : MO_NONE;
			else
				inst.ldst_order.store = (o0) ? MO_LO_RELEASE : MO_RELEASE;
			inst.ldst_order.rs = regRm(binst); // status register
		}

	common:
		inst.flags = set_mem_extend(inst.flags, size);
		if (size != SZ_X)
			inst.flags |= W32;

		inst.ldst.rt = regRd(binst);
		inst.ldst.rn = regRnSP(binst); // stack pointer can be base register
		break;
	}
	case Literal: {
		bool vector = (binst >> 26) & 1;

		inst.op = (vector) ? A64_LDR_FP : A64_LDR;
		if (vector) {
			switch (top2) { // opc
			case 0b00: inst.flags = set_mem_extend(inst.flags, FSZ_S); break;
			case 0b01: inst.flags = set_mem_extend(inst.flags, FSZ_D); break;
			case 0b10: inst.flags = set_mem_extend(inst.flags, FSZ_Q); break;
			case 0b11: return errinst("loads_and_stores/Literal: unallocated instruction");
			}
		} else {
			switch (top2) { // opc
			case 0b00: inst.flags = W32 | set_mem_extend(inst.flags, SZ_W); break;         // 32-bit LDR
			case 0b01: inst.flags = set_mem_extend(inst.flags, SZ_X); break;               // 64-bit LDR
			case 0b10: inst.flags = set_mem_extend(inst.flags, SXTW); break;               // LDRSW
			case 0b11: return errinst("loads_and_stores/Literal: PRFM not yet supported"); // XXX PRFM
			}
		}

		inst.offset = sext(((binst >> 5) & 0b1111111111111111111) << 2, 19+2); // imm19 * 4
		inst.ldst.rt = regRd(binst);
		break;
	}
	case PairNoAlloc: // fall through
	case Pair: {
		bool load = (binst >> 22) & 1;
		bool vector = (binst >> 26) & 1;
		int scale = 0; // The byte offset is stored as <imm7>/<scale>.

		if (vector) {
			inst.op = (load) ? A64_LDP_FP : A64_STP_FP;
			switch (top2) {
			case 0b00: scale =  4; inst.flags = set_mem_extend(inst.flags, FSZ_S); break; // Single (32 bits)
			case 0b01: scale =  8; inst.flags = set_mem_extend(inst.flags, FSZ_D); break; // Double (64 bits)
			case 0b10: scale = 16; inst.flags = set_mem_extend(inst.flags, FSZ_Q); break; // Quad  (128 bits)
			default:
				return errinst("loads_and_stores: invalid opc for LDP/STP");
			}
		} else {
			inst.op = (load) ? A64_LDP : A64_STP;
			switch (top2) {
			case 0b00: scale = 4; inst.flags = W32 | set_mem_extend(inst.flags, SZ_W); break; // 32-bit
			case 0b01: scale = 4; inst.flags = set_mem_extend(inst.flags, SXTW); break;       // LDPSW, STPSW
			case 0b10: scale = 8; inst.flags = set_mem_extend(inst.flags, SZ_X); break;       // 64-bit
			default:
				return errinst("loads_and_stores: invalid opc for LDP/STP");
			}
		}

		// PairNoAlloc is just Pair with a non-temporal hint, so no changes apart from the opcode.
		if (kind == PairNoAlloc) {
			switch (inst.op) {
			case A64_LDP:    inst.op = A64_LDNP; break;
			case A64_STP:    inst.op = A64_STNP; break;
			case A64_LDP_FP: inst.op = A64_LDNP_FP; break;
			case A64_STP_FP: inst.op = A64_STNP_FP; break;
			default:
				// Shut up impossible unhandled enum warning.
				return errinst("loads_and_stores/Pair: cannot happen");
			}
		}

		uint imm7 = (binst >> 15) & 0b1111111;
		inst.offset = scale * sext(imm7, 7);

		inst.ldst.rt = regRd(binst);
		inst.ldst.rn = regRnSP(binst); // stack pointer can be base register
		inst.ldst.rt2 = (binst >> 10) & 0b11111;
		break;
	}
	case UnscaledImm: // fall through -- like Register, but with different immediate (and LDUR/STUR opcodes)
	case Register: {
		bool vector = (binst >> 26) & 1;
		bool load = false;
		int scale = 0;    // for unsigned AM_OFF_IMM
		int shift = top2; // log2(accessed bits); shift amount for AM_OFF_REG and for AM_OFF_EXT

		switch (top2) {
		case 0b00: scale = 1; break;
		case 0b01: scale = 2; break;
		case 0b10: scale = 4; break;
		case 0b11: scale = 8; break;
		}

		u8 opc = (binst >> 22) & 0b11;
		if (vector) {
			FPSize size = 0;
			switch (opc) {
			case 0b00: load = false; size = top2; break; // STR_FP 8, 16, 32, 64 bit
			case 0b01: load = true;  size = top2; break; // LDR_FP ------
			case 0b10: // STR_FP 128 bit
				if (top2 != 0b00) {
					return errinst("loads_and_stores/Register: bad instruction");
				}
				load = false;
				scale = 16;
				shift = 4;
				size = FSZ_Q;
				break;
			case 0b11: // LDR_FP 128 bit
				if (top2 != 0b00) {
					return errinst("loads_and_stores/Register: bad instruction");
				}
				load = true;
				scale = 16;
				shift = 4;
				size = FSZ_Q;
				break;
			}

			if (kind == UnscaledImm) {
				inst.op = (load) ? A64_LDUR_FP : A64_STUR_FP;
			} else {
				inst.op = (load) ? A64_LDR_FP : A64_STR_FP;
			}
			inst.flags = set_mem_extend(inst.flags, size);
		} else {
			Size size = top2;
			bool sign_extend = false;
			bool w32 = (size != SZ_X);

			switch (opc) {
			case 0b00: load = false; sign_extend = 0; break; // store
			case 0b01: load = true;  sign_extend = 0; break; // load
			case 0b10: // load, sign-extend to 64 bit
				if (size == SZ_X) {
					// XXX PRFM -- prefetch memory
					return errinst("loads_and_stores/Register: PRFM not yet implemented");
				}
				load = true;
				sign_extend = 1;
				w32 = false;
				break;
			case 0b11: // load, sign-extend to 32 bit
				// 32-bit sign-extend to 32-bit does not make sense.
				// 64-bit sign-extend to 32-bit is similarly nonsensical.
				if (size >= SZ_W) {
					return errinst("loads_and_stores/Register: bad instruction");
				}
				load = true;
				sign_extend = 1;
				break;
			}

			if (kind == UnscaledImm) {
				inst.op = (load) ? A64_LDUR : A64_STUR;
			} else {
				inst.op = (load) ? A64_LDR : A64_STR;
			}
			if (w32) {
				inst.flags |= W32;
			}
			inst.flags = set_mem_extend(inst.flags, (sign_extend<<2) | size);
		}

		AddrMode mode = get_addrmode(inst.flags);
		switch (mode) {
		case AM_POST: // fall through
		case AM_PRE: {
			uint imm9 = (binst >> 12) & 0b111111111;
			inst.offset = sext(imm9, 9); // unscaled
			break;
		}
		case AM_OFF_IMM: {
			if (kind == UnscaledImm) {
				u32 imm9 = (binst >> 12) & 0b111111111;
				inst.offset = sext(imm9, 9);
				break;
			}

			// Normal LDR/STR scaled, unsigned immediate offset.
			uint imm12 = (binst >> 10) & 0b111111111111;
			inst.offset = scale * imm12;
			break;
		}
		case AM_OFF_REG: { // shifted register
			bool S = (binst >> 12) & 1;
			inst.shift.type = SH_LSL;
			inst.shift.amount = (S) ? shift : 0;
			inst.ldst.rm = regRm(binst);
			break;
		}
		case AM_OFF_EXT: { // extended register
			bool S = (binst >> 12) & 1;
			inst.extend.type = (binst >> 13) & 0b111; // option(3)
			inst.extend.lsl = (S) ? shift : 0;
			inst.ldst.rm = regRm(binst);
			break;
		}
		}

		inst.ldst.rt = regRd(binst);
		inst.ldst.rn = regRnSP(binst); // stack pointer can be base register
		break;
	}
	case Unprivileged:
		return errinst("loads_and_stores: Unprivileged instructions not supported");
	case Atomic: {
		u8 size = top2;
		u8 o3opc = (binst >> 12) & 0b1111; // o3 + opc

		bool acquire = (binst >> 23) & 1;
		bool release = (binst >> 22) & 1;

		// All < 64-bit accsizes store into the Wx registers; there are no
		// zero-extension or sign-extension options.
		if (size != SZ_X)
			inst.flags |= W32;

		switch (o3opc) {
		case 0b0000: inst.op = A64_LDADD;  break;
		case 0b0001: inst.op = A64_LDCLR;  break;
		case 0b0010: inst.op = A64_LDEOR;  break;
		case 0b0011: inst.op = A64_LDSET;  break;
		case 0b0100: inst.op = A64_LDSMAX; break;
		case 0b0101: inst.op = A64_LDSMIN; break;
		case 0b0110: inst.op = A64_LDUMAX; break;
		case 0b0111: inst.op = A64_LDUMIN; break;
		case 0b1000: inst.op = A64_SWP;    break;
		case 0b1100:
			// Load-AcquirePC Register -- ought to be in the Exclusive group,
			// but presumably did not fit in.
			inst.op = A64_LDR;
			inst.ldst_order.load = MO_ACQUIRE_PC;
			break;
		}

		inst.flags = set_mem_extend(inst.flags, size);
		inst.ldst_order.load = (acquire) ? MO_ACQUIRE : MO_NONE;
		inst.ldst_order.store = (release) ? MO_RELEASE : MO_NONE;
		inst.ldst.rt = regRd(binst);
		inst.ldst.rn = regRnSP(binst); // stack pointer can be base register
		inst.ldst.rs = regRm(binst);
		break;
	}
	case PAC:
		return errinst("loads_and_stores: Pointer Auth instructions not supported");
	}

	return inst;
}

// decode decodes n binary instructions from the input buffer and
// puts them in the output buffer, which must have space for n Insts.
int decode(u32 *in, uint n, Inst *out) {
	uint i;
	for (i = 0; i < n; i++) {
		u32 binst = in[i];
		u32 op0 = (binst >> 25) & 0b1111;
		switch (op0) {
		case 0b0000:
			out[i] = (Inst){
				op: A64_UDF,
			};
			break;
		case 0b1000:
		case 0b1001: // 100x
			out[i] = data_proc_imm(binst);
			break;
		case 0b1010:
		case 0b1011: // 101x
			out[i] = branches(binst);
			break;
		case 0b0100:
		case 0b0110:
		case 0b1100:
		case 0b1110: // x1x0
			out[i] = loads_and_stores(binst);
			break;
		case 0b0101:
		case 0b1101: // x101
			out[i] = data_proc_reg(binst);
			break;
		case 0b0111:
		case 0b1111: // x111
			out[i] = errinst("FP+SIMD processing not supported");
			break;
		default:
			out[i] = UNKNOWN_INST;
			break;
		}
	}
	return i;
}

enum {
	MAX_OBUF_SIZE = 1 * 1024 * 1024 * 1024 // 1 GB
};

// char *mnemonics[], generated using mnemonics.awk.
#include "mnemonics.h"

// main: read instructions as list of hexadecimal numbers, decode, print results.
int main(int argc, char **argv) {
	size_t capacity = 1;
	u64 ninst = 0;
	u32 *ibuf = calloc(capacity, sizeof(u32)); // dynamically resized
	Inst *obuf = NULL;

	char line[16]; // really, eight digits and some slack
	while (fgets(line, sizeof(line), stdin)) {
		ninst++;
		if (ninst >= capacity) {
			capacity *= 2;
			ibuf = realloc(ibuf, capacity * sizeof(u32));
			if (ibuf == NULL) {
				fprintf(stderr, "out of memory");
				return 1;
			}

			if (ninst * sizeof(Inst) >= MAX_OBUF_SIZE) {
				fprintf(stderr, "too many instructions, output buffer limit reached");
				return 2;
			}
		}

		sscanf(line, "%x", &ibuf[ninst-1]);
	}

	obuf = calloc(ninst, sizeof(Inst));
	if (obuf == NULL) {
		fprintf(stderr, "out of memory");
		return 1;
	}

	double ibufM = ((double)ninst*sizeof(u32)) / (1024.0 * 1024.0);
	double obufM = ((double)ninst*sizeof(Inst)) / (1024.0 * 1024.0);
	double totalM = ibufM + obufM;
	printf("#inst:     %7ld\nibuf size: %6.1fM\nobuf size: %6.1fM\ntotal:     %6.1fM\n",
		ninst, ibufM, obufM, totalM);

	decode(ibuf, ninst, obuf);

	for (uint i = 0; i < ninst; i++) {
		Inst inst = obuf[i];
		char regch = (inst.flags & W32) ? 'W' : 'X';
		char flagsch = (inst.flags & SET_FLAGS) ? 'S' : ' ';


		switch (inst.op) {
		case A64_UNKNOWN: printf("%04x ???\n", 4*i);                      continue;
		case A64_ERROR:   printf("%04x error \"%s\"\n", 4*i, inst.error); continue;
		case A64_UDF:     printf("%04x udf\n", 4*i);                       continue;
		default:
			break; // normal instruction
		}

		// We do not disambiguate here -- all instructions are printed
		// the same; for example, instructions with two immediates have
		// the imm field printed too.
		printf("%04x %-12s%c %c%d, %c%d, %c%d, imm=%lu, imm2=(%u,%u), offset=%+ld, w32=%o, set_flags=%o, memext=%o, addrmode=%o, cond=%x\n",
			4*i, mnemonics[inst.op], flagsch,
			regch, inst.rd, regch, inst.rn, regch, inst.rm,
			inst.imm, inst.bfm.lsb, inst.bfm.width, inst.offset, inst.flags&W32, inst.flags&SET_FLAGS,
			get_mem_extend(inst.flags), get_addrmode(inst.flags), get_cond(inst.flags));
	}

	return 0;
}
