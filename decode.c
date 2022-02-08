#include <stdbool.h>
#include <stdint.h>

#define FARMDEC_INTERNAL 1
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

Cond fad_get_cond(u8 flags) {
	return (flags >> 4) & 0b1111;
}

static u8 set_cond(u8 flags, Cond cond) {
	cond &= 0xF;
	flags &= 0x0F;
	return (cond << 4) | flags;
}

static u8 invert_cond(u8 flags) {
	Cond cond = fad_get_cond(flags);
	return set_cond(flags, cond ^ 0b001); // invert LSB
}

// Addressing mode, for Loads and Stores.
AddrMode fad_get_addrmode(u8 flags) {
	return (AddrMode)((flags >> 5) & 0b111);
}

static u8 set_addrmode(u8 flags, AddrMode mode) {
	return ((mode&0b111) << 5) | (flags&0b11111);
}

// How much memory to load/store (access size) and whether to sign-
// or zero-extend the value.
ExtendType fad_get_mem_extend(u8 flags) {
	return (ExtendType)((flags >> 2) & 0b111);
}

static u8 set_mem_extend(u8 flags, ExtendType memext) {
	return ((memext&0b111) << 2) | (flags&0b11100011);
}

VectorArrangement fad_get_vec_arrangement(u8 flags) {
	return (VectorArrangement)((flags >> 2) & 0b111);
}

static u8 set_vec_arrangement(u8 flags, VectorArrangement va) {
	return ((va&0b111) << 2) | (flags&0b11100011);
}

FPSize fad_get_prec(u8 flags) {
	return (FPSize)((flags >> 1) & 0b111);
}

static u8 set_prec(u8 flags, FPSize prec) {
	return ((prec&0b111) << 1) | (flags&0b11110001);
}

FPSize fad_size_from_vec_arrangement(VectorArrangement va) {
	return (FPSize)(va >> 1);
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
	s64 mask = (s64)1 << (b - 1);
	return (((s64)x) ^ mask) - mask;
}

// Helpers for data_proc_imm
static Inst find_bfm_alias(Op op, bool w32, Reg rd, Reg rn, u8 immr, u8 imms);
static u64 decode_bitmask(u8 immN, u8 imms, u8 immr, bool w32);

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
		u64 immhi = binst & (0b1111111111111111111 << 5);
		u64 immlo = top3 & 0b011;
		u64 uimm = (immhi >> (5-2)) | immlo; // pos(immhi) = 5; 2: len(immlo)

		u64 scale = (inst.op == A64_ADRP) ? 4096 : 1; // ADRP: Page (4K) Address
		s64 simm = scale * sext(uimm, 21);            // PC-relative → signed
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

		// ADDS/SUBS and thus CMN/CMP interpret R31 as the zero register,
		// while normal ADD and SUB treat it as the stack pointer.
		inst.rd = (inst.flags & SET_FLAGS) ? regRd(binst) : regRdSP(binst);
		inst.rn = regRnSP(binst);

		if (inst.rd == ZERO_REG && (inst.flags & SET_FLAGS)) {
			switch (inst.op) {
			case A64_ADD_IMM: inst.op = A64_CMN_IMM; break;
			case A64_SUB_IMM: inst.op = A64_CMP_IMM; break;
			default:
				break; // impossible; shuts up warning
			}
		} else if (inst.op == A64_ADD_IMM && !shift_by_12 && unshifted_imm == 0 &&
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

		u8 immr = (binst >> 16) & 0b111111;
		u8 imms = (binst >> 10) & 0b111111;
		u8 N = (inst.flags & W32) ? 0 : (binst >> 22) & 1; // N is part of imm for 64-bit variants
		inst.imm = decode_bitmask(N, imms, immr, inst.flags & W32);

		// ANDS and by extension TST interpret R31 as the zero register, while
		// regular immediate AND interprets it as the stack pointer.
		inst.rd = (inst.flags & SET_FLAGS) ? regRd(binst) : regRdSP(binst);
		inst.rn = regRn(binst);
		break;
	}
	case Move: {
		u8 hw = (binst >> 21) & 0b11;
		u8 shift = 16 * hw;
		u64 imm16 = (binst >> 5) & 0xFFFF;

		switch (top3 & 0b011) { // opc
		case 0b00: // MOVN: Move with NOT
			inst.op = A64_MOV_IMM;
			inst.imm = ~ (imm16 << shift);
			break;
		case 0b01:
			return UNKNOWN_INST;
		case 0b10: // MOVZ: zero other bits
			inst.op = A64_MOV_IMM;
			inst.imm = imm16 << shift;
			break;
		case 0b11: // MOVK: keep other bits
			inst.op = A64_MOVK;
			inst.movk.imm16 = imm16;
			inst.movk.lsl = shift;
			break;
		}

		inst.rd = regRd(binst);
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

static int fallback_highest_bit(unsigned int x) {
	// 43210  n
	// -----  -
	// 01001  0
	//  0100  1
	//   010  2
	//    01  3 <- our 0-based result
	int n = 0;
	do {
		x >>= 1;
		n++;
	} while (x != 1);
	return n;
}

// Returns the 0-based index of the highest bit. Should be compiled down
// to a single native instruction.
static int highest_bit(unsigned int x) {
#ifndef __has_builtin
	return fallback_highest_bit(x);
#else
#  if __has_builtin (__builtin_clz)
	// 31 - #leading zeroes
	// or 63 - #leading zeroes
	return (8*sizeof(x)-1) - __builtin_clz(x); // GCC and Clang.
#  else
	return fallback_highest_bit(x);
#  endif
#endif
}

// Rotate the len-bit number x n places to the right. Based on the first
// example at https://en.wikipedia.org/wiki/Bitwise_operation#Circular_shifts
// (except turned around, to make it rotare right).
static u64 ror(u64 x, uint n, uint len) {
	u64 raw = ((x >> n) | (x << (len - n)));
	if (len == 64)
		return raw;
	return raw & ((1uL<<len) - 1); // truncate left side to len bits
}

// Implementation of the A64 pseudocode function DecodeBitMasks (pp. 1683-1684).
//
// The logical immediate instructions encode 32-bit or 64-bit masks using merely
// 12 or 13 bits. We want the decoded mask in our Inst.imm field. We only need
// the "wmask" of DecodeBitMasks, so return only that.
static u64 decode_bitmask(u8 immN, u8 imms, u8 immr, bool w32) {
	uint M = (w32) ? 32 : 64;

	// Guarantee it's only the number of bits in the pseudocode signature.
	immN &= 1;
	imms &= 0b111111;
	immr &= 0b111111;

	// length of bitmask (1..6)
	uint len = highest_bit((immN << 6) | ((~imms)&0b111111));

	// 1..6 consecutive ones, basis of pattern
	uint levels = 0;
	for (uint i = 0; i < len; i++)
		levels = (levels << 1) | 1;

	uint S = imms & levels;
	uint R = immr & levels;
	uint esize = 1 << len;     // 2, 4, 8, 16, 32, 64

	// welem: pattern of 1s then zero-extended to esize
	// e.g. esize = 8; S+1 = 4 → welem = 0b00001111
	u64 welem = 0;
	for (uint i = 0; i < S+1; i++)
		welem = (welem << 1) | 1;

	// wmask = Replicate(ROR(welem, R));
	welem = ror(welem, R, esize);
	u64 wmask = 0;
	for (uint i = 0; i < M; i += esize) {
		wmask = (wmask << esize) | welem;
	}

	return wmask;
}

// There is no BFM, SBFM or UBFM instruction without a more specific alias.
// Returns a complete instruction, with flags and fields properly set.
static Inst find_bfm_alias(Op op, bool w32, Reg rd, Reg rn, u8 immr, u8 imms) {
	Inst inst = UNKNOWN_INST;
	u8 all_ones = (w32) ? 31 : 63; // for ASR, LSL, LSR
	u8 bits = (w32) ? 32 : 64;     // for BFC, BFI, LSL, SBFIZ, UBFIZ

	inst.rd = rd;
	inst.rn = rn;
	if (w32) {
		inst.flags |= W32;
	}

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
		// There is a very long op1 field where most of the values are the same
		// between the groups, forming a sort of prefix tree. Using the proper
		// order, we can differentiate the groups at the bits where the prefix
		// tree has branches.
		bool exc = !((binst >> 24) & 1); // exception → zero bit
		bool msr = (binst >> 20) & 1;    // system register move (MSR, MRS)
		bool sys = (binst >> 19) & 1;    // system instruction (SYS, SYSL)
		bool pstate = (binst >> 14) & 1; // PSTATE (MSR (imm),CFINV)

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
				return errinst("branches/Exceptions: bad instruction");
			}

			inst.imm = (binst >> 5) & 0xFFFF; // imm16
			break;
		}

		if (msr) {
			bool load = (binst >> 21) & 1;
			bool o0 = (binst >> 19) & 1;
			u8 o0val = (o0) ? 3 : 2;

			inst.op = (load) ? A64_MRS : A64_MSR_REG;
			inst.imm = (o0val << 14) | ((binst >> 5) & 0b11111111111111); // <o0val>:op1(3):CRn(4):CRm(4):op2(3)
			inst.rt = regRd(binst);
			break;
		}

		if (sys) {
			bool load = (binst >> 21) & 1;
			inst.op = (load) ? A64_SYSL : A64_SYS;
			inst.sys.op1 = (binst >> 16) & 0b111;
			inst.sys.op2 = (binst >> 5) & 0b111;
			inst.sys.crn = (binst >> 12) & 0b1111;
			inst.sys.crm = (binst >> 8) & 0b1111;
			inst.rt = regRd(binst);
			break;
		}

		if (pstate) {
			u8 op1 = (binst >> 16) & 0b111;
			u8 op2 = (binst >> 5) & 0b111;
			u8 crm = (binst >> 8) & 0b1111; // immediate

			switch ((op1 << 3) | op2) { // op1(3):op2(3)
			case 0b000000: inst.op = A64_CFINV;  return inst;
			case 0b000001: inst.op = A64_XAFlag; return inst;
			case 0b000010: inst.op = A64_AXFlag; return inst;
			case 0b000011: inst.msr_imm.psfld = PSF_UAO;     break;
			case 0b000100: inst.msr_imm.psfld = PSF_PAN;     break;
			case 0b000101: inst.msr_imm.psfld = PSF_SPSel;   break;
			case 0b011001: inst.msr_imm.psfld = PSF_SSBS;    break;
			case 0b011010: inst.msr_imm.psfld = PSF_DIT;     break;
			case 0b011110: inst.msr_imm.psfld = PSF_DAIFSet; break;
			case 0b011111: inst.msr_imm.psfld = PSF_DAIFClr; break;
			default:
				return errinst("branches/PSTATE: unknown MSR_IMM pstatefield");
			}
			inst.op = A64_MSR_IMM;
			inst.msr_imm.imm = crm;
			break;
		}

		// Longest prefix → Barriers or Hints, but we need to check the prefix!
		u16 op1 = (binst >> 12) & 0b11111111111111;
		if (op1 == 0b01000000110010) {
			inst.op = A64_HINT;
			// XXX decode further in the future
			return inst;
		}
		if (op1 != 0b01000000110011)
			return UNKNOWN_INST;

		u8 op2 = (binst >> 5) & 0b111;
		u8 crm = (binst >> 8) & 0b1111;
		switch (op2) {
		case 0b010: inst.op = A64_CLREX; inst.imm = crm; break;
		case 0b101: inst.op = A64_DMB; inst.imm = crm; break;
		case 0b110: inst.op = A64_ISB; inst.imm = crm; break;
		case 0b111:
			if (regRd(binst) != ZERO_REG)
				return errinst("branches/Barriers: bad instruction");
			inst.op = A64_SB;
			// CRm is always 0.
			break;
		case 0b100:
			switch (crm) {
			case 0b0000: inst.op = A64_SSBB; break;
			case 0b0100: inst.op = A64_PSSBB; break;
			default:
				inst.op = A64_DSB;
				inst.imm = crm;
				break;
			}
			break;
		default:
			return UNKNOWN_INST;
		}
		break;
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

		inst.offset = sext(((binst >> 5) & 0b1111111111111111111) << 2, 19+2); // imm19 * 4
		inst.rt = binst & 0b11111;
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
		inst.rt = binst & 0b11111;
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
		case 0b011: inst.op = A64_REV; break;
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

		// Unlike AddSubShifted, R31 is intepreted as the Stack Pointer,
		// unless it's ADDS/SUBS/CMN/CMP.
		inst.rd = (inst.flags & SET_FLAGS) ? regRd(binst) : regRdSP(binst);
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

		switch ((binst >> 10) & 0b111111) {
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
			inst.op = (immediate) ? A64_CCMN_IMM : A64_CCMN_REG;
		} else {
			inst.op = (immediate) ? A64_CCMP_IMM : A64_CCMP_REG;
		}

		inst.flags |= SET_FLAGS;
		inst.flags = set_cond(inst.flags, (binst >> 12) & 0b1111);
		inst.ccmp.nzcv = binst & 0b1111;

		inst.rn = regRn(binst);
		if (immediate) {
			inst.ccmp.imm5 = (binst >> 16) & 0b11111;
		} else {
			inst.rm = regRm(binst);
		}
		break;
	}
	case CondSelect: {
		// Combine the op bit (30) and the op2 field (11-10) to fully
		// determine the selection opcode.
		u32 op = (((binst >> 30) & 1) << 2) | ((binst >> 10) & 0b11);
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
		} else if (inst.rm == inst.rn && inst.op == A64_CSINV) {
			inst.op = A64_CINV;
			inst.flags = invert_cond(inst.flags);
		}
		break;
	}
	case DataProc3: {
		bool sub = (binst >> 15) & 1;
		switch ((binst >> 21) & 0b111) {
		case 0b000: inst.op = (sub) ? A64_MSUB : A64_MADD; break;
		case 0b001: inst.op = (sub) ? A64_SMSUBL : A64_SMADDL; break;
		case 0b010: inst.op = A64_SMULH; break;
		case 0b101: inst.op = (sub) ? A64_UMSUBL : A64_UMADDL; break;
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

		inst.flags = set_vec_arrangement(inst.flags, va);
		inst.simd_ldst.nreg = nreg;
		inst.rt = regRd(binst);
		inst.rn = regRnSP(binst);
		inst.rm = regRm(binst);

		// Post-indexed addrmode, either register-based or an immediate that
		// depends on the number of bytes read/written.
		if (fad_get_addrmode(inst.flags) == AM_POST) {
			if (inst.rm == ZERO_REG) {
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

		VectorArrangement va = (size << 1) | Q;

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

		inst.flags = set_vec_arrangement(inst.flags, va);
		inst.simd_ldst.nreg = nreg;
		inst.rt = regRd(binst);
		inst.rn = regRnSP(binst);

		// Post-indexed addrmode, either register-based or an immediate that
		// depends on the number of bytes read/written.
		if (fad_get_addrmode(inst.flags) == AM_POST) {
			inst.rm = regRm(binst);
			if (inst.rm == ZERO_REG) {
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
				inst.rs = regRm(binst);
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
				inst.rs = regRm(binst);
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
			inst.rt2 = (binst >> 10) & 0b11111;
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

		inst.rt = regRd(binst);
		inst.rn = regRnSP(binst); // stack pointer can be base register
		break;
	}
	case Literal: {
		bool vector = (binst >> 26) & 1;

		inst.op = (vector) ? A64_LDR_FP : A64_LDR;
		if (vector) {
			switch (top2) { // opc
			case 0b00: inst.flags = set_prec(inst.flags, FSZ_S); break;
			case 0b01: inst.flags = set_prec(inst.flags, FSZ_D); break;
			case 0b10: inst.flags = set_prec(inst.flags, FSZ_Q); break;
			case 0b11: return errinst("loads_and_stores/Literal: unallocated instruction");
			}
		} else {
			switch (top2) { // opc
			case 0b00: inst.flags = W32 | set_mem_extend(inst.flags, SZ_W); break; // 32-bit LDR
			case 0b01: inst.flags = set_mem_extend(inst.flags, SZ_X); break;       // 64-bit LDR
			case 0b10: inst.flags = set_mem_extend(inst.flags, SXTW); break;       // LDRSW
			case 0b11: inst.op = A64_PRFM; break;                                  // PRFM (literal)
			}
		}

		inst.offset = sext(((binst >> 5) & 0b1111111111111111111) << 2, 19+2); // imm19 * 4
		inst.rt = regRd(binst);
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
			case 0b00: scale =  4; inst.flags = set_prec(inst.flags, FSZ_S); break; // Single (32 bits)
			case 0b01: scale =  8; inst.flags = set_prec(inst.flags, FSZ_D); break; // Double (64 bits)
			case 0b10: scale = 16; inst.flags = set_prec(inst.flags, FSZ_Q); break; // Quad  (128 bits)
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

		inst.rt = regRd(binst);
		inst.rn = regRnSP(binst); // stack pointer can be base register
		inst.rt2 = (binst >> 10) & 0b11111;
		break;
	}
	case UnscaledImm: // fall through -- like Register, but with different immediate
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

			inst.op = (load) ? A64_LDR_FP : A64_STR_FP;
			inst.flags = set_prec(inst.flags, size);
		} else {
			Size size = top2;
			bool sign_extend = false;
			bool w32 = (size != SZ_X);

			switch (opc) {
			case 0b00: load = false; sign_extend = 0; break; // store
			case 0b01: load = true;  sign_extend = 0; break; // load
			case 0b10: // load, sign-extend to 64 bit
				if (size == SZ_X) { // PRFM (unsigned immediate, extended register or unscaled offset)
					inst.op = A64_PRFM;
					inst.rt = regRd(binst); // encodes prfop
					inst.rn = regRnSP(binst);

					if (fad_get_addrmode(inst.flags) == AM_OFF_EXT) {
						inst.extend.type = (binst >> 13) & 0b111;      // option(3)
						inst.extend.lsl = ((binst >> 12) & 1) ? 3 : 0; // S(1)
						inst.rm = regRmSP(binst);
					} else if (kind == UnscaledImm) {
						u32 imm9 = (binst >> 12) & 0b111111111;
						inst.offset = sext(imm9, 9);
					} else { // unsigned immediate
						uint imm12 = (binst >> 10) & 0b111111111111;
						inst.offset = 8 * imm12; // scale: 8
					}
					return inst;
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

			inst.op = (load) ? A64_LDR : A64_STR;
			if (w32) {
				inst.flags |= W32;
			}
			inst.flags = set_mem_extend(inst.flags, (sign_extend<<2) | size);
		}

		AddrMode mode = fad_get_addrmode(inst.flags);
		switch (mode) {
		default:
			return UNKNOWN_INST;

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
			inst.rm = regRm(binst);
			break;
		}
		case AM_OFF_EXT: { // extended register
			bool S = (binst >> 12) & 1;
			inst.extend.type = (binst >> 13) & 0b111; // option(3)
			inst.extend.lsl = (S) ? shift : 0;
			inst.rm = regRm(binst);
			break;
		}
		}

		inst.rt = regRd(binst);
		inst.rn = regRnSP(binst); // stack pointer can be base register
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
		inst.rt = regRd(binst);
		inst.rn = regRnSP(binst); // stack pointer can be base register
		inst.rs = regRm(binst);
		break;
	}
	case PAC:
		return errinst("loads_and_stores: Pointer Auth instructions not supported");
	}

	return inst;
}

static Inst scalar_floating_point(u32 binst);
static Inst decode_simd(u32 binst);

static Inst data_proc_float_and_simd(u32 binst) {
	Inst inst = UNKNOWN_INST;

	u8 op0 = (binst >> 28) & 0b1111;

	switch (op0) {
	case 0b0001:
	case 0b0011:
	case 0b1001:
	case 0b1011: // x0x1
		return scalar_floating_point(binst);
	case 0b1100:
		return errinst("cryptographic instructions not supported");
	case 0b1000:
	case 0b1010:
	case 0b1101:
	case 0b1110:
	case 0b1111:
		return errinst("unknown SIMD instruction");
	default:
		return decode_simd(binst);
	}

	return inst;
}

static double vfp_expand_imm(u8 imm);

// These instructions all use the vector registers as scalar floats of
// varying precision, and their encoding is rather similar to each other
// and rather different from the encoding of vector operations.
static Inst scalar_floating_point(u32 binst) {
	Inst inst = UNKNOWN_INST;

	bool b24 = (binst >> 24) & 1;
	bool b21 = (binst >> 21) & 1;
	u8 op3 = (binst >> 10) & 0b111111111;

	// Field "type" (bits 23-22): precision, but not in FPSize order
	switch ((binst >> 22) & 0b11) {
	case 0b00: inst.flags = set_prec(inst.flags, FSZ_S); break;
	case 0b01: inst.flags = set_prec(inst.flags, FSZ_D); break;
	case 0b10:
		// Do not give up just yet; FMOV_TOP2GPR and GRP2TOP have type=0b10.
		// In these cases, we overwrite inst with the right values again.
		inst = errinst("scalar_floating_point: bad precision");
		break;
	case 0b11: inst.flags = set_prec(inst.flags, FSZ_H); break;
	}

	enum {
		Unknown,     // b24 | b21 |    op3
		DataProc3,   //  1  |  x  | xxxxxxxxx
		ConvFixed,   //  0  |  0  | xxxxxxxxx
		ConvInt,     //  0  |  1  | xxx000000
		DataProc1,   //  0  |  1  | xxxx10000
		Compare,     //  0  |  1  | xxxxx1000
		Immediate,   //  0  |  1  | xxxxxx100
		CondCompare, //  0  |  1  | xxxxxxx01
		DataProc2,   //  0  |  1  | xxxxxxx10
		CondSelect,  //  0  |  1  | xxxxxxx11
	} kind = Unknown;

	if (b24) {
		kind = DataProc3;
	} else if (!b21) {
		kind = ConvFixed;
	} else if ((op3 & 0b111111) == 0b000000) {
		kind = ConvInt;
	} else if ((op3 & 0b11111)  ==  0b10000) {
		kind = DataProc1;
	} else if ((op3 & 0b1111)   ==   0b1000) {
		kind = Compare;
	} else if ((op3 & 0b111)    ==    0b100) {
		kind = Immediate;
	} else if ((op3 & 0b11)     ==     0b01) {
		kind = CondCompare;
	} else if ((op3 & 0b11)     ==     0b10) {
		kind = DataProc2;
	} else if ((op3 & 0b11)     ==     0b11) {
		kind = CondSelect;
	} else {
		kind = Unknown;
	}

	switch (kind) {
	case Unknown:
		return UNKNOWN_INST;
	case DataProc3: {
		bool o0 = (binst >> 15) & 1;
		bool o1 = (binst >> 21) & 1;

		switch ((o1 << 1) | o0) {
		case 0b00: inst.op = A64_FMADD;  break;
		case 0b01: inst.op = A64_FMSUB;  break;
		case 0b10: inst.op = A64_FNMADD; break;
		case 0b11: inst.op = A64_FNMSUB; break;
		}

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		inst.ra = (binst >> 10) & 0b11111;
		break;
	}
	case ConvFixed: {
		bool sf = (binst >> 31) & 1; // bit 31 is the size flag: 0→32, 1→64
		if (!sf) {
			inst.flags |= W32;
		}

		u8 opcode = (binst >> 16) & 0b111;
		switch (opcode) {
		case 0b000: inst.op = A64_FCVT_GPR; inst.fcvt.sgn = true; break;
		case 0b001: inst.op = A64_FCVT_GPR; inst.fcvt.sgn = false; break;
		case 0b010: inst.op = A64_CVTF; inst.fcvt.sgn = true; break;
		case 0b011: inst.op = A64_CVTF; inst.fcvt.sgn = false; break;
		default:
			errinst("scalar_floating_point/ConvFixed: bad opcode");
		}

		u8 scale = (binst >> 10) & 0b111111;

		inst.fcvt.mode = FPR_ZERO;
		inst.fcvt.fbits = 64 - scale; // always 64 - scale, even for w32
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;
	}
	case ConvInt: {
		bool sf = (binst >> 31) & 1; // bit 31 is the size flag: 0→32, 1→64
		if (!sf) {
			inst.flags |= W32;
		}

		u8 rmode = (binst >> 19) & 0b11;
		switch (rmode) {
		case 0b00: inst.fcvt.mode = FPR_TIE_EVEN; break;
		case 0b01: inst.fcvt.mode = FPR_POS_INF; break;
		case 0b10: inst.fcvt.mode = FPR_NEG_INF; break;
		case 0b11: inst.fcvt.mode = FPR_ZERO;    break;
		}

		u8 opcode = (binst >> 16) & 0b111;
		switch (opcode) {
		case 0b000: inst.op = A64_FCVT_GPR; inst.fcvt.sgn = true;  break;
		case 0b001: inst.op = A64_FCVT_GPR; inst.fcvt.sgn = false; break;
		case 0b010: inst.op = A64_CVTF; inst.fcvt.sgn = true;  break;
		case 0b011: inst.op = A64_CVTF; inst.fcvt.sgn = false; break;
		case 0b100: inst.op = A64_FCVT_GPR; inst.fcvt.sgn = true;  inst.fcvt.mode = FPR_TIE_AWAY; break; // FCVTAS
		case 0b101: inst.op = A64_FCVT_GPR; inst.fcvt.sgn = false; inst.fcvt.mode = FPR_TIE_AWAY; break; // FCVTAU
		case 0b110:
			// FMOV: Clear unused immediate again.
			switch (rmode) {
			case 0b00: inst.imm = 0; inst.op = A64_FMOV_VEC2GPR; break;
			case 0b01: inst.imm = 0; inst.op = A64_FMOV_TOP2GPR; break;
			case 0b11:
				inst.op = A64_FJCVTZS;
				inst.fcvt.mode = FPR_ZERO;
				inst.fcvt.sgn = true;
				break;
			default:
				return errinst("scalar_floating_point/ConvInt: bad rmode for opcode=0b110");
			}
			inst.rd = regRd(binst);
			inst.rn = regRn(binst);
			return inst;
		case 0b111:
			// FMOV: Clear unused immediate again.
			switch (rmode) {
			case 0b00: inst.imm = 0; inst.op = A64_FMOV_GPR2VEC; break;
			case 0b01: inst.imm = 0; inst.op = A64_FMOV_GPR2TOP; break;
			default:
				return errinst("scalar_floating_point/ConvInt: bad rmode for opcode=0b111");
			}
			inst.rd = regRd(binst);
			inst.rn = regRn(binst);
			return inst;
		}

		inst.fcvt.fbits = 0; // integer conversion → 0
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;
	}
	case DataProc1: {
		u8 opcode = (binst >> 15) & 0b111111;

		switch (opcode) {
		case 0b000000: inst.op = A64_FMOV_REG; break;
		case 0b000001: inst.op = A64_FABS; break;
		case 0b000010: inst.op = A64_FNEG; break;
		case 0b000011: inst.op = A64_FSQRT; break;
		case 0b000100: inst.op = A64_FCVT_S; break;
		case 0b000101: inst.op = A64_FCVT_D; break;
		case 0b000111: inst.op = A64_FCVT_H; break;
		case 0b001000: inst.op = A64_FRINT; inst.frint.mode = FPR_TIE_EVEN; break; // frintN
		case 0b001001: inst.op = A64_FRINT; inst.frint.mode = FPR_POS_INF;  break; // frintP
		case 0b001010: inst.op = A64_FRINT; inst.frint.mode = FPR_NEG_INF;  break; // frintM
		case 0b001011: inst.op = A64_FRINT; inst.frint.mode = FPR_ZERO;     break; // frintZ
		case 0b001100: inst.op = A64_FRINT; inst.frint.mode = FPR_TIE_AWAY; break; // frintA
			// 0b001101 is unallocated
		case 0b001110: inst.op = A64_FRINTX; inst.frint.mode = FPR_CURRENT; break; // frintX
		case 0b001111: inst.op = A64_FRINT;  inst.frint.mode = FPR_CURRENT; break; // frintI
		case 0b010000: inst.op = A64_FRINT;  inst.frint.mode = FPR_ZERO;    inst.frint.bits = 32; break; // frint32Z
		case 0b010001: inst.op = A64_FRINTX; inst.frint.mode = FPR_CURRENT; inst.frint.bits = 32; break; // frint32X
		case 0b010010: inst.op = A64_FRINT;  inst.frint.mode = FPR_ZERO;    inst.frint.bits = 64; break; // frint64Z
		case 0b010011: inst.op = A64_FRINTX; inst.frint.mode = FPR_CURRENT; inst.frint.bits = 64; break; // frint64X
		default:
			return errinst("scalar_floating_point/DataProc1: bad opcode");
		}

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;
	}
	case Compare: {
		u8 opcode = binst & 0b11111;
		switch (opcode) {
		case 0b00000: inst.op = A64_FCMP_REG;   break;
		case 0b01000: inst.op = A64_FCMP_ZERO;  break;
		case 0b10000: inst.op = A64_FCMPE_REG;  break;
		case 0b11000: inst.op = A64_FCMPE_ZERO; break;
		default:
			return errinst("scalar_floating_point/Compare: bad opcode");
		}

		inst.rn = regRn(binst);
		if (inst.op == A64_FCMP_REG || inst.op== A64_FCMPE_REG) {
			inst.rm = regRm(binst);
		}
		break;
	}
	case Immediate: {
		u8 imm8 = (binst >> 13) & 0xff;
		inst.op = A64_FMOV_IMM;
		inst.fimm = vfp_expand_imm(imm8);
		inst.rd = regRd(binst);
		break;
	}
	case CondCompare: {
		Cond cond = (binst >> 12) & 0b1111;
		bool signalling = (binst >> 4) & 1;

		inst.op = (signalling) ? A64_FCCMPE : A64_FCCMP;
		inst.flags = set_cond(inst.flags, cond);
		inst.ccmp.nzcv = binst & 0b1111;
		inst.ccmp.imm5 = 0;
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		break;
	}
	case DataProc2: {
		u8 opcode = (binst >> 12) & 0b1111;
		switch (opcode) {
		case 0b0000: inst.op = A64_FMUL;   break;
		case 0b0001: inst.op = A64_FDIV;   break;
		case 0b0010: inst.op = A64_FADD;   break;
		case 0b0011: inst.op = A64_FSUB;   break;
		case 0b0100: inst.op = A64_FMAX;   break;
		case 0b0101: inst.op = A64_FMIN;   break;
		case 0b0110: inst.op = A64_FMAXNM; break;
		case 0b0111: inst.op = A64_FMINNM; break;
		case 0b1000: inst.op = A64_FNMUL;  break;
		}
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		break;
	}
	case CondSelect: {
		Cond cond = (binst >> 12) & 0b1111;
		inst.op = A64_FCSEL;
		inst.flags = set_cond(inst.flags, cond);
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		break;
	}
	}

	return inst;
}

static Inst decode_simd(u32 binst) {
	Inst inst = UNKNOWN_INST;

	u8 size = (binst >> 22) & 0b11;
	bool scalar = (binst >> 28) & 1;
	bool U = (binst >> 29) & 1; // often – but not always – means Unsigned
	bool Q = (binst >> 30) & 1; // use full 128-bit vectors? (vector arrangement = size:Q)
	VectorArrangement va = (size << 1) | Q;

	// A sensible default that can be overridden by specific instructions where needed.
	inst.flags = set_vec_arrangement(inst.flags, va);

	// Although the tables in the manual have distinct groups for scalar instructions
	// (e.g. "Advanced SIMD scalar copy") and vector instructions (e.g. "Advanced SIMD
	// copy"), they share the same encoding "marker" and have the exact structure.
	// Bit 28 differentiates between the two: 0→vector, 1→scalar.
	enum {
		Unknown,
		// CryptoAES,
		// CryptoThreeRegSHA,
		// CryptoTwoRegSHA,
		// ScalarThreeSameFP16,
		// ScalarTwoRegMiscFP16,
		// ThreeSameFP16,
		// TwoRegMiscFP16,
		// CryptoThreeRegImm2,
		// CryptoThreeRegSHA512,
		// CryptoFourReg,
		// XAR,
		// CryptoTwoRegSHA512,
		TableLookup,
		Permute,
		Extract,
		Copy,           // includes ScalarCopy
		ThreeSameExtra, // includes ScalarThreeSameExtra
		TwoRegMisc,     // includes ScalarTwoRegMisc
		Reduce,         // ScalarPairwise, AcrossLanes
		ThreeDiff,      // includes ScalarThreeDiff
		ThreeSame,      // includes ScalarThreeSame
		ModifiedImm,
		ShiftByImm,     // includes ScalarShiftByImm
		IndexedElem,    // includes ScalarXIndexedElem
	} kind = Unknown;

	bool b10 = (binst >> 10) & 1;
	bool b11 = (binst >> 11) & 1;
	bool b15 = (binst >> 15) & 1;
	u8 op21 = (binst >> 21) & 0b1111; // <24:21> is a good starting point for differentiation
	u8 op17 = (binst >> 17) & 0b1111; // <20:17> is a secondary marker
	u8 immh = (binst >> 19) & 0b1111; // differentiates ShiftByImm and ModifiedImm

	switch (op21) {
	case 0b0000:
	case 0b0010:
	case 0b0100:
	case 0b0110: // 0xx0
		if (op21 == 0b0000 && U && !b15 && !b10)
			kind = Extract;
		else if (op21 == 0b0000 && !b15 && b10)
			kind = Copy;
		else if (op21 == 0b0000 && !b15 && !b11 && !b10)
			kind = TableLookup;
		else if (!b15 && b11 && !b10)
			kind = Permute;
		else if (b15 && b10)
			kind = ThreeSameExtra;
		break;
	case 0b0001:
	case 0b0011:
	case 0b0101:
	case 0b0111: // 0xx1
		if (op17 == 0b0000 && b11 && !b10)
			kind = TwoRegMisc;
		else if (op17 == 0b1000 && b11 && !b10)
			kind = Reduce;
		else if (!b11 && !b10)
			kind = ThreeDiff;
		else if (b10)
			kind = ThreeSame;
		break;
	case 0b1000:
	case 0b1001:
	case 0b1010:
	case 0b1011:
	case 0b1100:
	case 0b1101:
	case 0b1110:
	case 0b1111: // 1xxx
		if (immh == 0b0000 && b10)
			kind = ModifiedImm;
		else if (immh != 0b0000 && b10)
			kind = ShiftByImm;
		else if (!b10)
			kind = IndexedElem;
		break;
	}

	switch (kind) {
	case Unknown:
		return errinst("SIMD: bad instruction");
	case TableLookup: {
		bool tbx = (binst >> 12) & 1;
		u8 len = (binst >> 13) & 0b11; // #regs = len + 1
		inst.op = (tbx) ? A64_TBX : A64_TBL;
		inst.imm = len + 1;
		inst.rd = regRd(binst);
		inst.rn = regRn(binst); // table stored in Vn and successors (modulo 32)
		inst.rm = regRm(binst);
		break;
	}
	case Permute: {
		u8 opcode = (binst >> 12) & 0b111;
		switch (opcode) {
		case 0b001: inst.op = A64_UZP1; break;
		case 0b010: inst.op = A64_TRN1; break;
		case 0b011: inst.op = A64_ZIP1; break;
		case 0b101: inst.op = A64_UZP2; break;
		case 0b110: inst.op = A64_TRN2; break;
		case 0b111: inst.op = A64_ZIP2; break;
		default:
			return errinst("SIMD/Permute: bad opcode");
		}

		// Standard vector arrangement in size:Q, as set by the default handling
		// at the top of decode_simd.

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		break;
	}
	case Extract: {
		u8 imm4 = (binst >> 11) & 0b1111; // stores the index
		inst.op = A64_EXT;
		inst.flags = set_vec_arrangement(inst.flags, (Q) ? VA_16B : VA_8B);
		inst.imm = (Q) ? imm4 : (imm4 & 0b111);
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		break;
	}
	case Copy: {
		u8 imm4 = (binst >> 11) & 0b1111; // really not an immediate, more of an opcode
		u8 imm5 = (binst >> 16) & 0b11111;
		bool op = (binst >> 29) & 1;

		FPSize prec;
		u32 idx;
		if ((imm5 & 0b1111) == 0b1000) {      // x1000
			idx = imm5 >> 4;
			prec = FSZ_D;
		} else if ((imm5 & 0b111) == 0b100) { // xx100
			idx = imm5 >> 3;
			prec = FSZ_S;
		} else if ((imm5 & 0b11) == 0b10) {   // xxx10
			idx = imm5 >> 2;
			prec = FSZ_H;
		} else if ((imm5 & 0b1) == 1) {       // xxxx1
			idx = imm5 >> 1;
			prec = FSZ_B;
		} else {
			return errinst("SIMD/Copy: bad imm5 field / destination index");
		}

		inst.flags = set_vec_arrangement(inst.flags, (prec << 1) | Q);

		// INS_ELEM: uses imm4 as source index, not opcode.
		if (op) {
			inst.op = A64_INS_ELEM;
			inst.flags = set_vec_arrangement(inst.flags, (prec << 1) | 1); // always 128-bit
			inst.ins_elem.dst = idx;
			switch (prec) {
			case FSZ_B: inst.ins_elem.src = imm4 >> 0; break; // xxxx
			case FSZ_H: inst.ins_elem.src = imm4 >> 1; break; // xxx0
			case FSZ_S: inst.ins_elem.src = imm4 >> 2; break; // xx00
			case FSZ_D: inst.ins_elem.src = imm4 >> 3; break; // x000
			default:
				return errinst("SIMD/Copy: bad precision");
			}
			inst.rd = regRd(binst);
			inst.rn = regRn(binst);
			break;
		}

		// imm4 as opcode
		switch (imm4) {
		case 0b0000:
			inst.op = A64_DUP_ELEM;
			inst.flags |= (scalar) ? SIMD_SCALAR : 0;
			break;
		case 0b0001:
			inst.op = A64_DUP_GPR;
			inst.flags |= (prec != FSZ_D) ? W32 : 0;
			break;
		case 0b0011:
			inst.op = A64_INS_GPR;
			inst.flags = set_vec_arrangement(inst.flags, (prec << 1) | 1); // always 128-bit
			inst.flags |= (prec != FSZ_D) ? W32 : 0;
			break;
		case 0b0101: {
			inst.op = A64_SMOV;
			bool is_long = ((imm5 & 0b10000) > 0);
			inst.flags = set_vec_arrangement(inst.flags, (prec << 1) | is_long);
			inst.flags |= (!Q) ? W32 : 0;
			inst.flags |= SIMD_SIGNED;
			break;
		}
		case 0b0111: {
			inst.op = A64_UMOV;
			bool is_long = ((imm5 & 0b10000) > 0);
			inst.flags = set_vec_arrangement(inst.flags, (prec << 1) | is_long);
			inst.flags |= (!Q) ? W32 : 0;
			break;
		}
		default:
			return errinst("SIMD/Copy: bad instruction");
		}

		inst.imm = idx;
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;
	}
	case ThreeSameExtra: {
		u8 opcode = (binst >> 11) & 0b1111;
		switch (opcode) {
		case 0b0000: inst.op = A64_SQRDMLAH_VEC; inst.flags |= (scalar) ? SIMD_SCALAR : 0; break;
		case 0b0001: inst.op = A64_SQRDMLSH_VEC; inst.flags |= (scalar) ? SIMD_SCALAR : 0; break;
		case 0b0010: inst.op = A64_DOT_VEC;      inst.flags |= (U) ? 0 : SIMD_SIGNED;      break;
		case 0b1000:
		case 0b1001:
		case 0b1010:
		case 0b1011: // 10xx (xx = rot)
			inst.op = A64_FCMLA_VEC;
			switch (opcode & 0b11) {
			case 0b00: inst.imm =   0; break;
			case 0b01: inst.imm =  90; break;
			case 0b10: inst.imm = 180; break;
			case 0b11: inst.imm = 270; break;
			}
			break;
		case 0b1100:
		case 0b1110: // 11x0 (x = rot)
			inst.op = A64_FCADD;
			inst.imm = ((opcode >> 1) & 1) ? 270 : 90;
			break;
		}

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		break;
	}
	case TwoRegMisc: {
		u8 opcode = (binst >> 12) & 0b11111;
		bool set_signedness = false; // set Inst.flags.signed = NOT (<U bit>)?
		bool set_scalarity = false;  // set Inst.flags.scalar = <scalar bit>?

		bool is_fp = false;

		switch (opcode) {
		case 0b00000: inst.op = (U) ? A64_REV32_VEC : A64_REV64_VEC; break;
		case 0b00001: inst.op = A64_REV16_VEC; break;
		case 0b00010: inst.op = A64_ADDLP; set_signedness = 1; break;
		case 0b00011: inst.op = (U) ? A64_USQADD : A64_SUQADD; set_scalarity = 1; break;
		case 0b00100: inst.op = (U) ? A64_CLZ_VEC : A64_CLS_VEC; break;
		case 0b00101:
			inst.op = (U) ? ((size==0) ? A64_NOT_VEC : A64_RBIT_VEC) : A64_CNT;
			inst.flags = set_vec_arrangement(inst.flags, (FSZ_B << 1) | Q);
			break;
		case 0b00110: inst.op = A64_ADALP; set_signedness = 1; break;
		case 0b00111: inst.op = (U) ? A64_SQNEG : A64_SQABS; inst.flags |= SIMD_SIGNED; set_scalarity = 1; break;
		case 0b01000: inst.op = (U) ? A64_CMGE_ZERO : A64_CMGT_ZERO; set_scalarity = 1; break;
		case 0b01001: inst.op = (U) ? A64_CMLE_ZERO : A64_CMEQ_ZERO; set_scalarity = 1; break;
		case 0b01010: inst.op = A64_CMLT_ZERO; set_scalarity = 1; break;
		case 0b01011: inst.op = (U) ? A64_NEG_VEC : A64_ABS_VEC; set_scalarity = 1; break;
		case 0b01100: is_fp = true; inst.op = (U) ? A64_FCMGE_ZERO : A64_FCMGT_ZERO; set_scalarity = 1; break;
		case 0b01101: is_fp = true; inst.op = (U) ? A64_FCMLE_ZERO : A64_FCMEQ_ZERO; set_scalarity = 1; break;
		case 0b01110: is_fp = true; inst.op = A64_FCMLT_ZERO; set_scalarity = 1; break;
		case 0b01111: is_fp = true; inst.op = (U) ? A64_FNEG_VEC : A64_FABS_VEC; set_scalarity = 1; break;
		case 0b10010: inst.op = (U) ? A64_SQXTUN : A64_XTN; set_scalarity = 1; break;
		case 0b10100: inst.op = A64_QXTN; set_signedness = 1; set_scalarity = 1; break;
		case 0b10110:
			inst.op = (U) ? A64_FCVTXN : A64_FCVTN;
			set_scalarity = 1;

			// We store vector arrangement of the narrower destination, like all other
			// Narrowing instructions. Since we set the arrangement here, is_fp = false.
			switch ((binst >> 22) & 1) {
			case 0: inst.flags = set_vec_arrangement(inst.flags, (FSZ_H << 1) | Q); break;
			case 1: inst.flags = set_vec_arrangement(inst.flags, (FSZ_S << 1) | Q); break;
			}

			break;
		case 0b10111:
			// FCVTL{2} decodes the source vector's size in size(0):Q a bit differently compared
			// to most FP operations, so we set the arrangement manually and have is_fp = false.
			inst.op = A64_FCVTL;
			switch ((size&1) << 1 | Q) {
			case 0b00: inst.flags = set_vec_arrangement(inst.flags, VA_4H); break;
			case 0b01: inst.flags = set_vec_arrangement(inst.flags, VA_8H); break;
			case 0b10: inst.flags = set_vec_arrangement(inst.flags, VA_2S); break;
			case 0b11: inst.flags = set_vec_arrangement(inst.flags, VA_4S); break;
			}
			break;
		case 0b11000:
			is_fp = true;
			inst.op = A64_FRINT_VEC;
			inst.frint.mode = (U) ? FPR_TIE_AWAY : ((size & 0b10) ? FPR_POS_INF : FPR_TIE_EVEN);
			break;
		case 0b11001:
			is_fp = true;
			if (U) {
				inst.op = (size & 0b10) ? A64_FRINT_VEC : A64_FRINTX_VEC;
				inst.frint.mode = FPR_CURRENT;
			} else {
				inst.op = A64_FRINT_VEC;
				inst.frint.mode = (size & 0b10) ? FPR_ZERO : FPR_NEG_INF;
			}
			break;
		case 0b11010: is_fp = true; inst.op = A64_FCVT_VEC; inst.fcvt.mode = (size & 0b10) ? FPR_POS_INF : FPR_TIE_EVEN; set_signedness = 1; set_scalarity = 1; break;
		case 0b11011: is_fp = true; inst.op = A64_FCVT_VEC; inst.fcvt.mode = (size & 0b10) ? FPR_ZERO : FPR_NEG_INF; set_signedness = 1; set_scalarity = 1; break;
		case 0b11100:
			is_fp = true; // URECPE, URSQRTE not actually FP, but still only use the size<0> bit to encode size
			if (size & 0b10) {
				inst.op = (U) ? A64_URSQRTE : A64_URECPE;
			} else {
				inst.op = A64_FCVT_VEC;
				inst.fcvt.mode = FPR_TIE_AWAY;
				set_signedness = 1;
				set_scalarity = 1;
			}
			break;
		case 0b11101:
			is_fp = true;
			if (size & 0b10) {
				inst.op = (U) ? A64_FRSQRTE : A64_FRECPE;
			} else {
				inst.op = A64_CVTF_VEC;
				set_signedness = 1;
				set_scalarity = 1;
			}
			break;
		case 0b11110:
			is_fp = true;
			inst.op = (U) ? A64_FRINTX_VEC : A64_FRINT_VEC;
			inst.frint.bits = 32;
			inst.frint.mode = (U) ? FPR_CURRENT : FPR_ZERO;
			break;
		case 0b11111:
			is_fp = true;
			if (size & 0b10) {
				inst.op = (U) ? A64_FSQRT_VEC : A64_FRECPX;
				set_scalarity = 1;
			} else {
				inst.op = (U) ? A64_FRINTX_VEC : A64_FRINT_VEC;
				inst.frint.bits = 64;
				inst.frint.mode = (U) ? FPR_CURRENT : FPR_ZERO;
			}
			break;
		default:
			return errinst("SIMD/TwoRegMisc: bad opcode");
		}

		if (set_signedness) {
			inst.flags |= (U) ? 0 : SIMD_SIGNED;
			if (inst.op == A64_FCVT_VEC || inst.op == A64_CVTF_VEC)
				inst.fcvt.sgn = !U;
		}

		if (set_scalarity)
			inst.flags |= (scalar) ? SIMD_SCALAR : 0;

		// The floating-point instructions use the upper bit of size(2) to differentiate
		// instructions, leaving only the lower bit to toggle between single and double
		// precision.
		if (is_fp)
			inst.flags = set_vec_arrangement(inst.flags, (((size&1) ? FSZ_D : FSZ_S) << 1) | Q);

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;
	}
	case Reduce: {
		u8 opcode = (binst >> 12) & 0b11111;
		bool set_signedness = false; // set Inst.flags.signed = NOT (<U bit>)?
		VectorArrangement va = -1;
		bool is_fmin = (size & 0b10) != 0;

		inst.flags |= (scalar) ? SIMD_SCALAR : 0;

		switch (opcode) {
		case 0b00011: inst.op = A64_ADDLV; set_signedness = 1; break;
		case 0b01010: inst.op = A64_MAXV;  set_signedness = 1; break;
		case 0b11010: inst.op = A64_MINV;  set_signedness = 1; break;
		case 0b11011: inst.op = (scalar) ? A64_ADDP : A64_ADDV; break;
		case 0b01100:
			if (scalar) {
				inst.op = (is_fmin) ? A64_FMINNMP : A64_FMAXNMP;
				if (U) {
					va = (size&1) ? VA_2D : VA_2S;
				} else {
					va = VA_4H; // actually 2H (which does not exist)
				}
			} else {
				inst.op = (is_fmin) ? A64_FMINNMV : A64_FMAXNMV;
				va = (((U) ? ((size&1) ? FSZ_D : FSZ_S) : FSZ_H) << 1) | Q;
			}

			break;
		case 0b01101:
			inst.op = A64_FADDP;
			if (U) {
				va = (size&1) ? VA_2D : VA_2S;
			} else {
				va = VA_4H; // actually 2H (which does not exist)
			}
			break;
		case 0b01111:
			if (scalar) {
				inst.op = (is_fmin) ? A64_FMINP : A64_FMAXP;
				if (U) {
					va = (size&1) ? VA_2D : VA_2S;
				} else {
					va = VA_4H; // actually 2H (which does not exist)
				}
			} else {
				inst.op = (is_fmin) ? A64_FMINV : A64_FMAXV;
				va = (((U) ? ((size&1) ? FSZ_D : FSZ_S) : FSZ_H) << 1) | Q;
			}
			break;
		}

		if (set_signedness)
			inst.flags |= (U) ? 0 : SIMD_SIGNED;

		if (va != -1)
			inst.flags = set_vec_arrangement(inst.flags, va);

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;
	}
	case ThreeDiff: {
		u8 opcode = (binst >> 12) & 0b1111;
		bool set_signedness = false; // set Inst.flags.signed = NOT (<U bit>)?
		bool set_scalarity = false;  // set Inst.flags.scalar = <scalar bit>?

		switch (opcode) {
		case 0b0000: inst.op = A64_ADDL; set_signedness = true; break;
		case 0b0001: inst.op = A64_ADDW; set_signedness = true; break;
		case 0b0010: inst.op = A64_SUBL; set_signedness = true; break;
		case 0b0011: inst.op = A64_SUBW; set_signedness = true; break;
		case 0b0100: inst.op = A64_ADDHN; inst.flags |= (U) ? SIMD_ROUND : 0; break;
		case 0b0101: inst.op = A64_ABAL; set_signedness = true; break;
		case 0b0110: inst.op = A64_SUBHN; inst.flags |= (U) ? SIMD_ROUND : 0; break;
		case 0b0111: inst.op = A64_ABDL; set_signedness = true; break;
		case 0b1000: inst.op = A64_MLAL_VEC; set_signedness = true; break;
		case 0b1001: inst.op = A64_SQDMLAL_VEC; inst.flags |= SIMD_SIGNED; set_scalarity = true; break;
		case 0b1010: inst.op = A64_MLSL_VEC; set_signedness = true; break;
		case 0b1011: inst.op = A64_SQDMLAL_VEC; inst.flags |= SIMD_SIGNED; set_scalarity = true; break;
		case 0b1100: inst.op = A64_MULL_VEC; set_signedness = true; break;
		case 0b1101: inst.op = A64_SQDMULL_VEC; inst.flags |= SIMD_SIGNED; set_scalarity = true; break;
		case 0b1110: inst.op = A64_PMULL; break;
		default:
			return errinst("SIMD/ThreeDiff: bad obcode");
		}

		if (set_signedness)
			inst.flags |= (U) ? 0 : SIMD_SIGNED;

		if (set_scalarity)
			inst.flags |= (scalar) ? SIMD_SCALAR : 0;

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst); 
		break;
	}
	case ThreeSame: {
		u8 opcode = (binst >> 11) & 0b11111;
		bool set_signedness = false; // set Inst.flags.signed = NOT (<U bit>)?
		bool set_scalarity = false;  // set Inst.flags.scalar = <scalar bit>?

		bool is_fp = false;

		switch (opcode) {
		case 0b00000: inst.op = A64_HADD; set_signedness = 1; break;
		case 0b00001: inst.op = A64_QADD; set_signedness = 1; set_scalarity = 1; break;
		case 0b00010: inst.op = A64_HADD; set_signedness = 1; inst.flags |= SIMD_ROUND; break;
		case 0b00011:
			// size acts as sub-opcode
			switch (size) {
			case 0b00: inst.op = (U) ? A64_EOR_VEC : A64_AND_VEC;     break;
			case 0b01: inst.op = (U) ? A64_BSL     : A64_BIC_VEC_REG; break;
			case 0b10: inst.op = (U) ? A64_BIT     : A64_ORR_VEC_REG; break;
			case 0b11: inst.op = (U) ? A64_BIF     : A64_ORN_VEC;     break;
			}
			inst.flags = set_vec_arrangement(inst.flags, (FSZ_B << 1) | Q);
			break;
		case 0b00100: inst.op = A64_HSUB; set_signedness = 1; break;
		case 0b00101: inst.op = A64_QSUB; set_signedness = 1; set_scalarity = 1; break;
		case 0b00110: inst.op = (U) ? A64_CMHI_REG : A64_CMGT_REG; set_scalarity = 1; break;
 		case 0b00111: inst.op = (U) ? A64_CMHS_REG : A64_CMGE_REG; set_scalarity = 1; break;
		case 0b01000: inst.op = A64_SHL_REG; set_signedness = 1; set_scalarity = 1; break;
		case 0b01001: inst.op = A64_QSHL_REG; set_signedness = 1; set_scalarity = 1; break;
		case 0b01010: inst.op = A64_SHL_REG; set_signedness = 1; set_scalarity = 1; inst.flags |= SIMD_ROUND; break;
		case 0b01011: inst.op = A64_QSHL_REG; set_signedness = 1; set_scalarity = 1; inst.flags |= SIMD_ROUND; break;
		case 0b01100: inst.op = A64_MAX_VEC; set_signedness = 1; break;
		case 0b01101: inst.op = A64_MIN_VEC; set_signedness = 1; break;
		case 0b01110: inst.op = A64_ABD; set_signedness = 1; break;
		case 0b01111: inst.op = A64_ABA; set_signedness = 1; break;
		case 0b10000: inst.op = (U) ? A64_SUB_VEC : A64_ADD_VEC; set_scalarity = 1; break;
		case 0b10001: inst.op = (U) ? A64_CMEQ_REG : A64_CMTST; set_scalarity = 1; break;
		case 0b10010: inst.op = (U) ? A64_MLS_VEC : A64_MLA_VEC; break;
		case 0b10011: inst.op = (U) ? A64_PMUL : A64_MUL_VEC; break;
		case 0b10100: inst.op = A64_MAXP; set_signedness = 1; break;
		case 0b10101: inst.op = A64_MINP; set_signedness = 1; break;
		case 0b10110: inst.op = A64_SQDMULH_VEC; set_scalarity = 1; inst.flags |= (U) ? SIMD_ROUND : 0; break;
		case 0b10111: inst.op = A64_ADDP_VEC; break;
		case 0b11000: is_fp = true; inst.op = (U) ? ((size&0b10) ? A64_FMINNMP_VEC : A64_FMAXNMP_VEC) : ((size&0b10) ? A64_FMINNM_VEC : A64_FMAXNM_VEC); break;
		case 0b11001: is_fp = true; inst.op = (U) ? ((size&0b10) ? A64_FMLSL2_VEC : A64_FMLAL2_VEC) : ((size&0b10) ? A64_FMLS_VEC : A64_FMLA_VEC); break;
		case 0b11010:
			is_fp = true;
			inst.op = (U) ? ((size&0b10) ? A64_FABD_VEC : A64_FADDP_VEC) : ((size&0b10) ? A64_FSUB_VEC : A64_FADD_VEC);
			if (inst.op == A64_FABD_VEC) {
				set_scalarity = 1;
			}
			break;
		case 0b11011: is_fp = true; inst.op = (U) ? A64_FMUL_VEC : ((scalar) ? A64_FMULX : A64_FMULX_VEC); break;
		case 0b11100: is_fp = true; inst.op = (U) ? ((size&0b10) ? A64_FCMGT_REG : A64_FCMGE_REG) : A64_FCMEQ_REG; set_scalarity = 1; break;
		case 0b11101: is_fp = true; inst.op = (U) ? ((size&0b10) ? A64_FACGT : A64_FACGE) : ((size&0b10) ? A64_FMLSL_VEC : A64_FMLAL_VEC); break; set_scalarity = 1; break;
		case 0b11110: is_fp = true; inst.op = (U) ? ((size&0b10) ? A64_FMINP_VEC : A64_FMAXP_VEC) : ((size&0b10) ? A64_FMIN_VEC : A64_FMAX_VEC); break;
		case 0b11111: is_fp = true; inst.op = (U) ? A64_FDIV_VEC : ((size&0b10) ? A64_FRSQRTS_VEC : A64_FRECPS_VEC); break;
		}

		if (set_signedness)
			inst.flags |= (U) ? 0 : SIMD_SIGNED;

		if (set_scalarity)
			inst.flags |= (scalar) ? SIMD_SCALAR : 0;

		// The floating-point instructions use the upper bit of size(2) to differentiate
		// instructions, leaving only the lower bit to toggle between single and double
		// precision.
		if (is_fp)
			inst.flags = set_vec_arrangement(inst.flags, (((size&1) ? FSZ_D : FSZ_S) << 1) | Q);

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		break;
	}
	case ModifiedImm: {
		bool op = (binst >> 29) & 1;
		u8 cmode = (binst >> 12) & 0b1111;

		u8 abc = (binst >> 16) & 0b111;
		u8 defgh = (binst >> 5) & 0b11111;
		u64 imm8 = (abc << 5) | defgh; // a:b:c:d:e:f:g:h

		uint zlsl = 0; // zero-shifting left shift amount
		switch ((cmode >> 1) & 0b11) { // cmode<2:1>
		case 0b00: zlsl =  0; break;
		case 0b01: zlsl =  8; break;
		case 0b10: zlsl = 16; break;
		case 0b11: zlsl = 24; break;
		}

		uint olsl = (cmode & 1) ? 16 : 8; // ones-shifting left shift amount

		switch (cmode) {
		case 0b0000:
		case 0b0010:
		case 0b0100:
		case 0b0110: // 0xx0 → MOVI(32), MVNI(22)
			inst.op = A64_MOVI;
			inst.flags = set_vec_arrangement(inst.flags, (FSZ_S << 1) | Q);
			inst.imm = imm8 << zlsl;
			inst.imm = (op) ? ~inst.imm : inst.imm; // op=1 → invert
			break;
		case 0b1000:
		case 0b1010: // 10x0 → MOVI(16), MVNI(16)
			inst.op = A64_MOVI;
			inst.flags = set_vec_arrangement(inst.flags, (FSZ_H << 1) | Q);
			inst.imm = imm8 << zlsl;
			inst.imm = (op) ? ~inst.imm : inst.imm; // op=1 → invert
			break;
		case 0b1100:
		case 0b1101: // 110x → MOVI(32, shifting ones), MVNI(32, shifting ones)
			inst.op = A64_MOVI;
			inst.flags = set_vec_arrangement(inst.flags, (FSZ_S << 1) | Q);
			inst.imm = (imm8 << olsl) | ~(~((u64)0) << olsl); // shift #olsl ones
			inst.imm = (op) ? ~inst.imm : inst.imm; // op=1 → invert
			break;
		case 0b1110: // 1110 → MOVI(8), MOVI(64, vec), MOVI(64, scalar)
			inst.op = A64_MOVI;
			if (op == 0) {
				inst.flags = set_vec_arrangement(inst.flags, (FSZ_B << 1) | Q);
				inst.imm = imm8;
			} else {
				inst.flags |= (Q) ? 0 : SIMD_SCALAR;
				inst.flags = set_vec_arrangement(inst.flags, (FSZ_D << 1) | Q);

				// 8 x a : 8 x b : ... : 8 x h
				u64 abyte = (imm8 & 0b10000000) ? 0xff : 0;
				u64 bbyte = (imm8 & 0b01000000) ? 0xff : 0;
				u64 cbyte = (imm8 & 0b00100000) ? 0xff : 0;
				u64 dbyte = (imm8 & 0b00010000) ? 0xff : 0;
				u64 ebyte = (imm8 & 0b00001000) ? 0xff : 0;
				u64 fbyte = (imm8 & 0b00000100) ? 0xff : 0;
				u64 gbyte = (imm8 & 0b00000010) ? 0xff : 0;
				u64 hbyte = (imm8 & 0b00000001) ? 0xff : 0;
				inst.imm = (abyte << 56) | (bbyte << 48) | (cbyte << 40) | (dbyte << 32)
				         | (ebyte << 24) | (fbyte << 16) | (gbyte <<  8) | (hbyte <<  0);
			}
			break;
		case 0b0001:
		case 0b0011:
		case 0b0101:
		case 0b0111: // 0xx1 → ORR(32), BIC(32)
			inst.op = (op) ? A64_BIC_VEC_IMM : A64_ORR_VEC_IMM;
			inst.flags = set_vec_arrangement(inst.flags, (FSZ_S << 1) | Q);
			inst.imm = imm8 << zlsl;
			break;
		case 0b1001:
		case 0b1011: // 10x1 → ORR(16), BIC(16)
			inst.op = (op) ? A64_BIC_VEC_IMM : A64_ORR_VEC_IMM;
			inst.flags = set_vec_arrangement(inst.flags, (FSZ_H << 1) | Q);
			inst.imm = imm8 << zlsl;
			break;
		case 0b1111: {
			bool o2 = (binst >> 11) & 1;

			inst.op = A64_FMOV_VEC;
			if (!op && !o2)
				inst.flags = set_vec_arrangement(inst.flags, (FSZ_S << 1) | Q);
			else if (!op && o2)
				inst.flags = set_vec_arrangement(inst.flags, (FSZ_H << 1) | Q);
			else if (op && !o2)
				inst.flags = set_vec_arrangement(inst.flags, (FSZ_D << 1) | Q);
			else
				return errinst("SIMD/ModifiedImm: FMOV_VEC: unknown (op, o2) combination");

			inst.fimm = vfp_expand_imm(imm8);
			break;
		}
		default:
			return errinst("SIMD/ModifiedImm: bad cmode");
		}

		inst.rd = regRd(binst);
		break;
	}
	case ShiftByImm: {
		u8 opcode = (binst >> 11) & 0b11111;
		bool set_signedness = true; // set Inst.flags.signed = NOT (<U bit>)?
		bool set_scalarity = true;  // set Inst.flags.scalar = <scalar bit>?
		bool left_shift = false;

		switch (opcode) {
		case 0b00000: inst.op = A64_SHR; break;
		case 0b00010: inst.op = A64_SRA; break;
		case 0b00100: inst.op = A64_SHR; inst.flags |= SIMD_ROUND; break;
		case 0b00110: inst.op = A64_SRA; inst.flags |= SIMD_ROUND; break;
		case 0b01000: inst.op = A64_SRI; break;
		case 0b01010: inst.op = (U) ? A64_SLI : A64_SHL_IMM; left_shift = 1; set_signedness = 0; break;
		case 0b01110: inst.op = A64_QSHL_IMM; left_shift = 1; break;
		case 0b10000:
		case 0b10001: // 1000x (x = round?)
			inst.op = (U) ? A64_SQSHRUN : A64_SHRN;
			set_signedness = 0;
			set_scalarity = 0;
			inst.flags |= (opcode & 1) ? SIMD_ROUND : 0;
			break;
		case 0b10010: inst.op = A64_QSHRN; break;
		case 0b10011: inst.op = A64_QSHRN; inst.flags |= SIMD_ROUND; break;
		case 0b10100: inst.op = A64_SHLL; left_shift = 1; set_scalarity = 0; break;
		case 0b11100: inst.op = A64_CVTF_VEC; break; // fixed-point
		case 0b11111: inst.op = A64_FCVT_VEC; break; // fixed-point
		default:
			return errinst("SIMD/ShiftByImm: bad opcode");
		}

		if (set_signedness)
			inst.flags |= (U) ? 0 : SIMD_SIGNED;

		if (set_scalarity)
			inst.flags |= (scalar) ? SIMD_SCALAR : 0;

		u8 immh = (binst >> 19) & 0b1111;   // immh -- encodes size
		u8 immh_immb = (binst >> 16) & 0b1111111; // immh:immb -- encodes immediate (shift amount)
		FPSize size;

		if (immh == 0b0000)
			return errinst("SIMD/ShiftByImm: internal decoding to kind failed: immh=0000 → ModifiedImm");
		else if ((immh & 0b1111) == 0b0001)
			size = FSZ_B;
		else if ((immh & 0b1110) == 0b0010)
			size = FSZ_H;
		else if ((immh & 0b1100) == 0b0100)
			size = FSZ_S;
		else if ((immh & 0b1000) == 0b1000)
			size = FSZ_D;
		else
			return errinst("SIMD/ShiftByImm: bad immediate");

		uint imm = 0; // shift amount (or fractional bits for fixed-point CVTF, FCVT)
		switch (size) {
		case FSZ_B: imm = (left_shift) ? immh_immb -  8 :  16 - immh_immb; break;
		case FSZ_H: imm = (left_shift) ? immh_immb - 16 :  32 - immh_immb; break;
		case FSZ_S: imm = (left_shift) ? immh_immb - 32 :  64 - immh_immb; break;
		case FSZ_D: imm = (left_shift) ? immh_immb - 64 : 128 - immh_immb; break;
		default:
			return errinst("SIMD/ShiftByImm: bad immediate");
		}

		if (inst.op == A64_CVTF_VEC || inst.op == A64_FCVT_VEC) {
			inst.fcvt.mode = FPR_ZERO;
			inst.fcvt.fbits = imm;
			inst.fcvt.sgn = !U;
		} else {
			inst.imm = imm;
		}

		inst.flags = set_vec_arrangement(inst.flags, (size << 1) | Q);
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;
	}
	case IndexedElem: {
		u8 opcode = (binst >> 12) & 0b1111;
		bool set_signedness = false; // set Inst.flags.signed = NOT (<U bit>)?
		bool set_scalarity = false;  // set Inst.flags.scalar = <scalar bit>?

		// Bits H:L:M encode the index. Note that H is at pos 11, while
		// the lower two bits L:M are at pos 21:20.
		u8 hlm = ((binst >> (11-2)) & 0b100) | ((binst >> 20) & 0b11); 

		// use only 4 bits to encode Rm, restricting it to V0..V15?
		// used when the M bit is part of the index and not part of the register number
		bool short_rm = false;

		// The index usually depends on the size, with the exception of some instructions.
		u8 idx = 0;
		switch (size) {
		case 0b01: idx = hlm; short_rm = true; break; // FSZ_H -> H:L:M [0,1,2,3,4,5,6,7]
		case 0b10: idx = hlm >> 1; break;             // FSZ_S -> H:L   [0,1,2,3]
		case 0b11: idx = hlm >> 2; break;             // FSZ_D -> H     [0,1]
		default:
			break;
		}

		inst.imm = idx; // sensible default

		// FCMLA is an outlier: it has two immediates, the index and the rotation,
		// the latter of which is encoded in the middle of the opcode field!
		if (U && (opcode & 0b1001) == 0b0001) {
			inst.op = A64_FCMLA_ELEM;
			switch ((opcode>>1) & 0b11) {
			case 0b00: inst.fcmla_elem.rot =   0; break;
			case 0b01: inst.fcmla_elem.rot =  90; break;
			case 0b10: inst.fcmla_elem.rot = 180; break;
			case 0b11: inst.fcmla_elem.rot = 270; break;
			}

			// We can only index the lower halves.
			switch (size) {
			case 0b01: inst.fcmla_elem.idx = hlm >> 1; break; // H:L
			case 0b10: inst.fcmla_elem.idx = hlm >> 2; break; // H
			}

			inst.rd = regRd(binst);
			inst.rn = regRn(binst);
			inst.rm = regRm(binst);
			break;
		}

		switch (opcode) {
		case 0b0000:
			if (U) {
				inst.op = A64_MLA_ELEM;
				break;
			}

			inst.op = A64_FMLAL_ELEM;
			inst.flags = set_vec_arrangement(inst.flags, FSZ_H << 1 | Q);
			inst.imm = hlm;
			short_rm = true; // M is used in index
			break;
		case 0b0001: inst.op = A64_FMLA_ELEM; set_scalarity  = true; break;
		case 0b0010: inst.op = A64_MLAL_ELEM; set_signedness = true; break;
		case 0b0011: inst.op = A64_SQDMLAL_ELEM; set_scalarity = true; inst.flags |= SIMD_SIGNED; break;
		case 0b0100:
			if (U) {
				inst.op = A64_MLS_ELEM;
				break;
			}

			inst.op = A64_FMLSL_ELEM;
			inst.flags = set_vec_arrangement(inst.flags, FSZ_H << 1 | Q);
			inst.imm = hlm;
			short_rm = true; // M is used in index
			break;
		case 0b0101: inst.op = A64_FMLS_ELEM; set_scalarity  = true; break;
		case 0b0110: inst.op = A64_MLSL_ELEM; set_signedness = true; break;
		case 0b0111: inst.op = A64_SQDMLSL_ELEM; set_scalarity = true; inst.flags |= SIMD_SIGNED; break;
		case 0b1000:
			if (!U) {
				inst.op = A64_MUL_ELEM;
				break;
			}

			inst.op = A64_FMLAL2_ELEM;
			inst.flags = set_vec_arrangement(inst.flags, FSZ_H << 1 | Q);
			inst.imm = hlm;
			short_rm = true; // M is used in index
			break;
		case 0b1001: inst.op = (U) ? A64_FMULX_ELEM : A64_FMUL_ELEM; set_scalarity = true; break;
		case 0b1010: inst.op = A64_MULL_ELEM; set_signedness = true; break;
		case 0b1011: inst.op = A64_SQDMULL_ELEM; set_scalarity = true; inst.flags |= SIMD_SIGNED; break;
		case 0b1100:
			if (!U) {
				inst.op = A64_SQDMULH_ELEM;
				set_scalarity = true;
				inst.flags |= SIMD_SIGNED;
				break;
			}
			inst.op = A64_FMLSL2_ELEM;
			inst.flags = set_vec_arrangement(inst.flags, FSZ_H << 1 | Q);
			inst.imm = hlm;
			short_rm = true; // M is used in index
			break;
		case 0b1101:
			inst.op = (U) ? A64_SQRDMLAH_ELEM : A64_SQDMULH_ELEM; // SQRDMULH, SQRDMLAH
			inst.flags |= SIMD_SIGNED | SIMD_ROUND;
			set_scalarity = true;
			break;
		case 0b1110: // e.g.: udot v0.2s, v1.8b, v2.4b[idx] (idx = 0..3)
			inst.op = A64_DOT_ELEM;
			inst.flags = set_vec_arrangement(inst.flags, FSZ_B << 1 | Q);
			inst.imm = hlm >> 1; // H:L
			set_signedness = true;
			break;
		case 0b1111: inst.op = A64_SQRDMLSH_ELEM; inst.flags |= SIMD_SIGNED | SIMD_ROUND; set_scalarity = true; break;
		}

		if (set_signedness)
			inst.flags |= (U) ? 0 : SIMD_SIGNED;

		if (set_scalarity)
			inst.flags |= (scalar) ? SIMD_SCALAR : 0;

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		inst.rm = regRm(binst);
		inst.rm = (short_rm) ? (inst.rm&0b1111) : inst.rm;
		break;
	}
	}

	return inst;
}

// Expands an FMOV 8-bit FP immediate. Implementation of the pseudocode
//
//     bits(N) VFPExpandImm(bits(8) imm8)
//
// found in the ARM A64 ISA manual, specialised to bits(64). The imm8
// is constructed as follows:
//
//     7 6 5 4 3 2 1 0
//     ± -exp- --frc--
//
static double vfp_expand_imm(u8 imm8) {
	u64 sign = (((u64)imm8) >> 7) & 1;

	// imm8(exp) → "sign extend" to double exponent of 11 bits.
	u64 exp = (((u64)imm8) >> 4) & 0b111;
	u64 expsgn = (exp >> 2) & 1;
	u64 inverted_expsgn = (~expsgn) & 1; 
	u64 dexp = (inverted_expsgn<<10) | (expsgn<<9) | (expsgn<<8) | (expsgn<<7) | (expsgn<<6)
	            | (expsgn<<5) | (expsgn<<4) | (expsgn<<3) | (expsgn<<2) | (exp&0b11);

	// imm8(frac) → top 4 bits of double mantissa of 52 bits.
	u64 dfrac = (((u64)imm8) & 0b1111) << (52 - 4);

	union {
		u64    ieee754;
		double d;
	} foo;
	foo.ieee754 = (sign << 63) | ((dexp&0b11111111111) << 52) | dfrac;
	return foo.d;
}

// decode decodes n binary instructions from the input buffer and
// puts them in the output buffer, which must have space for n Insts.
int fad_decode(u32 *in, uint n, Inst *out) {
	Inst inst;
	uint i;
	for (i = 0; i < n; i++) {
		u32 binst = in[i];
		u32 op0 = (binst >> 25) & 0b1111;
		switch (op0) {
		case 0b0000:
			inst = (Inst){
				op: A64_UDF,
			};
			break;
		case 0b1000:
		case 0b1001: // 100x
			inst = data_proc_imm(binst);
			break;
		case 0b1010:
		case 0b1011: // 101x
			inst = branches(binst);
			break;
		case 0b0100:
		case 0b0110:
		case 0b1100:
		case 0b1110: // x1x0
			inst = loads_and_stores(binst);
			break;
		case 0b0101:
		case 0b1101: // x101
			inst = data_proc_reg(binst);
			break;
		case 0b0111:
		case 0b1111: // x111
			inst = data_proc_float_and_simd(binst);
			break;
		default:
			inst = UNKNOWN_INST;
			break;
		}

		if (inst.op == A64_UNKNOWN) {
			inst.imm = binst;
		}
		out[i] = inst;
	}
	return i;
}
