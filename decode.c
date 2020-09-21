#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

typedef unsigned int uint;
typedef uint8_t   u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t  s64;

typedef enum Op Op;
typedef u8 Reg;
typedef struct Inst Inst;

// Opcodes ordered and grouped according to the Top-level Encodings
// of the A64 Instruction Set Architecture (ARMv8-A profile) document,
// pages 1406-1473.
//
// Unlike normal ARM assembly, instructions are not overloaded here.
// Instead, there are markers like _IMM (immediate) and _SH (shifted
// register) where disambiguation is needed. Aliases are listed after
// the base instruction.
//
// However, whether the instruction sets flags or uses 32-bit registers
// is encoded in fields of Inst, not in the opcode, so ADDS is represented
// by A64_ADD_*.
enum Op {
	A64_UNKNOWN, // an invalid or (to us) unknown instruction
	A64_UDF,     // throws undefined exception

	/*** Data Processing -- Immediate ***/

	// PC-rel. addressing
	A64_ADR,     // ADR Xd, label  -- Xd ← PC + label
	A64_ADRP,    // ADRP Xd, label -- Xd ← PC + (label * 4K)

	// Add/subtract (immediate, with tags)
	A64_ADDG,
	A64_SUBG,

	// Add/subtract (immediate)
	A64_ADD_IMM,
	A64_SUB_IMM,

	// Logical (immediate)
	A64_AND_IMM,
	A64_ORR_IMM,
	A64_EOR_IMM,
	A64_TST_IMM, // TST Rd -- ANDS alias (Rd := RZR); normal ANDS is represented by A64_AND_IMM

	// XXX imm. moves

	// Bitfield
	A64_SBFM,
	A64_ASR_IMM,
	A64_SBFIZ,
	A64_SBFX,
	A64_SXTB,
	A64_SXTH,
	A64_SXTW,
	A64_BFM,
	A64_BFC,
	A64_BFI,
	A64_BFXIL,
	A64_UBFM,
	A64_LSL_IMM,
	A64_LSR_IMM,
	A64_UBFIZ,
	A64_UBFX,
	A64_UXTB,
	A64_UXTH,

	// Extract
	A64_EXTR,
	A64_ROR_IMM, // ROR Rd, Rs, #shift -- EXTR alias (Rm := Rs, Rn := Rs)

	/*** Branches, Exception Generating and System Instructions ***/

	A64_BCOND, // XXX idk, how to represent conds? all as different opcodes?

	// XXX system instructions

	// Unconditional branch (register)
	A64_BR,
	A64_BLR,
	A64_RET,

	// Unconditional branch (immediate)
	A64_B,
	A64_BL,

	// Compare and branch (immediate)
	A64_CBZ,
	A64_CBNZ,

	// Test and branch (immediate)
	A64_TBZ,
	A64_TBNZ
};

// XXX keep it at 16 bytes if possible
struct Inst {
	Op  op;
	u8 flags; // lower four bits: see enum flagmasks; upper four bit: see enum Cond
	Reg rd;  // destination register - Rd
	Reg rn;  // first (or only) operand - Rn, Rt
	Reg rm;  // second operand - Rm
	union {
		u64 imm;     // primary immediate
		s64 offset;  // alternative to imm for branches: PC-relative byte offset
		Reg ra;      // third operand for 3-source data proc instrs (MADD, etc.)
		u32 imm2[2]; // two immediates where necessary
	};
};

// The meaning of the Inst.imm2 field.
enum imm2_indexes {
	IMMR = 0,
	IMMS = 1
};

enum flagmasks {
	W32 = 1 << 0,       // use the 32-bit W0...W31 facets?
	SET_FLAGS = 1 << 1, // modify the NZCV flags? (S mnemonic suffix)
};

typedef enum Cond Cond;

// The condition bits used by conditial branches, selects and compares, stored in the
// upper four bit of the Inst.flags field. The first three bits determine the condition
// proper while the LSB inverts the condition if set.
enum Cond {
	COND_EQ = 0b0000,  // =
	COND_NE = 0b0001,  // ≠
	COND_CS = 0b0010,  // Carry Set
	COND_HS = COND_CS, // ≥, Unsigned
	COND_CC = 0b0011,  // Carry Clear
	COND_LO = COND_CC, // <, Unsigned
	COND_MI = 0b0100,  // < 0 (MInus)
	COND_PL = 0b0101,  // ≥ 0 (PLus)
	COND_VS = 0b0110,  // Signed Overflow
	COND_VC = 0b0111,  // No Signed Overflow
	COND_HI = 0b1000,  // >, Unsigned
	COND_LS = 0b1001,  // ≤, Unsigned
	COND_GE = 0b1010,  // ≥, Signed
	COND_LT = 0b1011,  // <, Signed
	COND_GT = 0b1100,  // >, Signed
	COND_LE = 0b1101,  // ≤, Signed
	COND_AL = 0b1110,  // Always true
	COND_NV = 0b1111,  // Always true (not "never" as in A32!)
};

static Cond get_cond(u8 flags) {
	return (flags >> 4) & 0b1111;
}

static u8 set_cond(u8 flags, Cond cond) {
	flags &= 0x0F; // clear cond bits first
	flags |= cond << 4;
	return flags;
}

static Inst UNKNOWN_INST = {
	op: A64_UNKNOWN,
	// all other fields: zero
};

// The destination register Rd, if present, occupies bits 0..4.
static Reg regRd(u32 binst) {
	return binst & 0b11111;
}

// The first operand register Rn, if present, occupies bits 5..9.
static Reg regRn(u32 binst) {
	return (binst >> 5) & 0b11111;
}

enum {
	// Register W31/X31 is a zero register for some instructions and the
	// stack pointer for others.
	ZERO_REG = 31
};

// sext sign-extends the b-bits number in x to 64 bit. The upper (64-b) bits
// must be zero. Seldom needed, but fiddly.
//
// Taken from https://graphics.stanford.edu/~seander/bithacks.html#VariableSignExtend
static s64 sext(u64 x, u8 b) {
	s64 mask = 1U << (b - 1);
	return (((s64)x) ^ mask) - mask;
}

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

		inst.flags &= ~W32; // no 32-bit variant for these

		u64 immhi = binst & (0b111111111111111111 << 5);
		u64 immlo = top3 & 0b011;
		inst.imm = (immhi >> (5-2)) | immlo; // pos(immhi) = 5; 2: len(immlo)

		inst.rd = regRd(binst);
		break;
	}
	case AddSubTags:
		if (top3 == 0b100)
			inst.op = A64_ADDG;
		else if (top3 == 0b110)
			inst.op = A64_SUBG;
		else
			return UNKNOWN_INST;

		printf("ADDG, SUBG not supported"); // XXX
		return UNKNOWN_INST;

	case AddSub: {
		bool is_add = (top3 & 0b010) == 0;
		inst.op = (is_add) ? A64_ADD_IMM : A64_SUB_IMM;
		if (top3 & 0b001)
			inst.flags |= SET_FLAGS;

		u64 unshifted_imm = (binst >> 10) & 0b111111111111;
		bool shift_by_12 = (binst & (1 << 22)) > 0;
		inst.imm = (shift_by_12) ? unshifted_imm << 12 : unshifted_imm;

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
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

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;
	}
	case Move: {
		printf("Move (imm) not supported");
		break; // XXX
	}
	case Bitfield: {
		switch (top3 & 0b011) { // base instructions at this point
		case 0b00: inst.op = A64_SBFM; break;
		case 0b01: inst.op = A64_BFM;  break;
		case 0b10: inst.op = A64_UBFM; break;
		default:
			return UNKNOWN_INST;
		}

		u64 immr = (binst >> 16) & 0b111111;
		u64 imms = (binst >> 10) & 0b111111;

		// Now that we have the immediates, we can check for fitting aliases.
		switch (inst.op) {
		case A64_SBFM:
			if (((inst.flags & W32) && imms == 31) || (!(inst.flags & W32) && imms == 63)) {
				// SBFM Rd, Rn, #immr, #31/#63 -> ASR (shift := immr)
				inst.op = A64_ASR_IMM;
				inst.imm = immr;
				goto regs;
			}
			if (imms < immr) {
				inst.op = A64_SBFIZ; // SBFIZ Rd, Rn, #lsb, #width
				inst.imm2[IMMR] = 0xDEAD;  // XXX immr = -lsb MOD 32/64 → lsb = ?
				inst.imm2[IMMS] = imms + 1; // imms = width - 1 → width = imms + 1
				goto regs;
			}
			if (immr == 0) {
				switch (imms) {
				case  7: inst.op = A64_SXTB; inst.imm = 0; goto regs; // no immediate
				case 15: inst.op = A64_SXTH; inst.imm = 0; goto regs; // no immediate
				case 31: inst.op = A64_SXTW; inst.imm = 0; goto regs; // no immediate
				}
			}
			inst.op = A64_SBFX; // SBFX Rd, Rn, #lsb, #width
			inst.imm2[IMMR] = immr;            // immr = lsb
			inst.imm2[IMMS] = imms - immr + 1; // imms = lsb + width - 1 → width = imms - lsb + 1;

			break;
			// There is no SBFM general case not handled by an alias
			// XXX is that true?
		case A64_BFM:
			printf("BFM not implemented"); // XXX
			break;

		case A64_UBFM:
			printf("UBFM not implemented"); // XXX
			// XXX how to re-use parts of SBFM?
			break;
		default:
			// cannot happen
			return UNKNOWN_INST;
		}

	regs:
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break; // XXX
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

	// For branches, the actual byte offset is the immediate times 4.
	// Hence the multiplications.

	// XXX virtual address would help in pre-calculating absolute address
	// XXX like fadec. We could do that in an after-pass to avoid an extra
	// XXX argument for every function.

	switch (kind) {
	case Unknown:
		return UNKNOWN_INST;

	case CondBranch:
		inst.flags = set_cond(inst.flags, binst & 0b1111);
		inst.offset = 4 * sext((binst >> 5) & 0b1111111111111111111, 19); // imm19
		break;

	case System:
		printf("System instructions not supported\n");
		break; // XXX _some_ should be decoded. barriers, CFINV

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
		inst.offset = 4 * sext(binst & 0b11111111111111111111111111, 26); // imm26
		break;

	case CmpAndBranch: {
		if ((top7 & 0b1000000) == 0)
			inst.flags |= W32;

		bool zero = (binst & (1 << 24)) == 0;
		inst.op = (zero) ? A64_CBZ : A64_CBNZ;

		inst.offset = 4 * sext((binst >> 5) & 0b1111111111111111111, 19); // imm19
		inst.rn = binst & 0b11111; // Rt; not modified → Inst.rn, not .rd
		break;
	}
	case TestAndBranch: {
		if ((top7 & 0b1000000) == 0)
			inst.flags |= W32;

		bool zero = (binst & (1 << 24)) == 0;
		inst.op = (zero) ? A64_TBZ : A64_TBNZ;

		inst.imm2[0] = 4 * sext((binst >> 5) & 0b11111111111111, 14); // imm14
		u32 b40 = (binst >> 19) & 0b11111;
		u32 b5 = binst & (1<<31);
		inst.imm2[1] = (b5 >> (31-5)) | b40; // b5:b40
		inst.rn = binst & 0b11111;           // Rt; not modified → Inst.rn, not .rd
		break;
	}
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
			printf("0x%04x: Loads & stores not supported\n", 4*i); // XXX
			out[i] = UNKNOWN_INST;
			break;
		case 0b0101:
		case 0b1101: // x101
			printf("0x%04x: Register data processing not supported\n", 4*i); // XXX
			out[i] = UNKNOWN_INST;
			break;
		case 0b0111:
		case 0b1111: // x111
			printf("0x%04x: FP+SIMD processing not supported\n", 4*i); // XXX
			out[i] = UNKNOWN_INST;
			break;
		default:
			out[i] = UNKNOWN_INST;
			break;
		}
	}
	return i;
}

// Buffers not in main because allocating hundreds of MB on the stack
// leads to a segfault.
#define NINST (28)
u32 ibuf[NINST];
Inst obuf[NINST];

int main(int argc, char **argv) {
	double ibufM = ((double)sizeof(ibuf)) / (1024.0 * 1024.0);
	double obufM = ((double)sizeof(obuf)) / (1024.0 * 1024.0);
	double totalM = ibufM + obufM;
	printf("#inst:     %6dM\nibuf size: %6.1fM\nobuf size: %6.1fM\ntotal:     %6.1fM\n",
		NINST / (1024 * 1024), ibufM, obufM, totalM);

	/*
	.LBB1_9:
		0000 and     w11, w20, #0xff
		0004 orr     x8, x8, x22
		0008 mov     x1, x19
		000c ldp     x20, x19, [sp, #48]
		0010 ldp     x22, x21, [sp, #32]
		0014 ldp     x24, x23, [sp, #16]
		0018 orr     x8, x8, x9
		001c orr     x8, x8, x10
		0020 orr     x0, x8, x11, lsl #48
		0024 ldp     x29, x30, [sp], #64
		0028 ret
	.LBB1_10:
		002c tst     w21, #0x40000000
		0030 mov     w8, #6
		0034 ubfx    w10, w21, #10, #12
		0038 mvn     w9, w21
		003c cinc    x22, x8, ne
		0040 lsl     x8, x10, #12
		0044 tst     w21, #0x400000
		0048 mov     w0, w21
		004c ubfx    w23, w21, #29, #1
		0050 lsr     w24, w9, #31
		0054 csel    x19, x10, x8, eq
		0058 mov     w20, w0
		005c mov     w0, w21
		0060 lsl     x10, x0, #56
		0064 lsl     x9, x23, #40
		0068 lsl     x8, x24, #32
		006c b       .LBB1_9

	*/
	unsigned char sample[] = {
		0x8b, 0x1e, 0x00, 0x12,
		0x08, 0x01, 0x16, 0xaa,
		0xe1, 0x03, 0x13, 0xaa,
		0xf4, 0x4f, 0x43, 0xa9,
		0xf6, 0x57, 0x42, 0xa9,
		0xf8, 0x5f, 0x41, 0xa9,
		0x08, 0x01, 0x09, 0xaa,
		0x08, 0x01, 0x0a, 0xaa,
		0x00, 0xc1, 0x0b, 0xaa,
		0xfd, 0x7b, 0xc4, 0xa8,
		0xc0, 0x03, 0x5f, 0xd6,

		0xbf, 0x02, 0x02, 0x72,
		0xc8, 0x00, 0x80, 0x52,
		0xaa, 0x56, 0x0a, 0x53,
		0xe9, 0x03, 0x35, 0x2a,
		0x16, 0x05, 0x88, 0x9a,
		0x48, 0xcd, 0x74, 0xd3,
		0xbf, 0x02, 0x0a, 0x72,
		0xe0, 0x03, 0x15, 0x2a,
		0xb7, 0x76, 0x1d, 0x53,
		0x38, 0x7d, 0x1f, 0x53,
		0x53, 0x01, 0x88, 0x9a,
		0xf4, 0x03, 0x00, 0x2a,
		0xe0, 0x03, 0x15, 0x2a,
		0x0a, 0x1c, 0x48, 0xd3,
		0xe9, 0x5e, 0x58, 0xd3,
		0x08, 0x7f, 0x60, 0xd3,
		0xe5, 0xff, 0xff, 0x17
	};

	// We just repeat the sample one instruction at a time until NINST is reached.
	unsigned char *p = (unsigned char *)sample;
	for (uint i = 0; i < NINST; i++) {
		// Little endian to native endianness. The values must be unsigned bytes.
		u32 binst = (p[0] << 0) | (p[1] << 8) | (p[2] << 16) | (p[3] << 24);
		ibuf[i] = binst;

		p += 4;
		if (p[0] == 0 && p[1] == 0 && p[2] == 0 && p[3] == 0) // end of sample array
			p = (unsigned char *) sample; // -> start from the beginning
	}

	decode(ibuf, NINST, obuf);

	if (1) {
		for (uint i = 0; i < NINST; i++) {
			Inst inst = obuf[i];
			char regch = (inst.flags & W32) ? 'W' : 'X';
			char flagsch = (inst.flags & SET_FLAGS) ? 'S' : ' ';

			// We do not disambiguate here -- all instructions are printed
			// the same; for example, instructions with two immediates (using
			// the imm2[] field) have the imm field printed too.
			printf("0x%04x: %d%c %c%d, %c%d, %c%d, imm=%lu, imm2=(%u,%u), offset=%4ld\n", 4*i, inst.op, flagsch,
				regch, inst.rd, regch, inst.rn, regch, inst.rm, inst.imm, inst.imm2[0], inst.imm2[1], inst.offset);
		}
	}

	return 0;
}
