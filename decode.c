#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;

typedef enum Op Op;
typedef u8 Reg;
typedef struct Inst Inst;

// Opcodes ordered and grouped according to the Top-level Encodings
// of the A64 Instruction Set Architecture (ARMv8-A profile) document,
// pages 1406-1473.
//
// Unlike normal ARM assembly, instructions are not overloaded here.
// Instead, there are markers like _IMM (immediate) and _SH (shifted
// register) where disambiguation is needed.
//
// However, whether the instruction sets flags or uses 32-bit registers
// is encoded in fields of Inst, not in the opcode, so ADDS is represented
// by A64_ADD_*.
enum Op {
	A64_UNKNOWN, // an invalid or (to us) unknown instruction
	A64_UDF,     // throws undefined exception

	A64_ADR,     // ADR Xd, label  -- Xd ← PC + label
	A64_ADRP,    // ADRP Xd, label -- Xd ← PC + (label * 4K)

	A64_ADDG,
	A64_SUBG,
	A64_ADD_IMM,
	A64_SUB_IMM,
};

struct Inst {
	Op   op;
	bool w32;       // use 32-bit register facet?
	bool set_flags; // set NZCV flags?
	Reg  rd;
	Reg  rn;
	u64  imm; // primary immediate
};

static Inst UNKNOWN_INST = {
	op: A64_UNKNOWN,
	// all other fields: zero
};

// The destination register Rd, if present, occupies bits 0..4.
static Reg regRd(u32 inst) {
	return inst & 0b11111;
}

// The source register Rn, if present, occupies bits 5..9.
static Reg regRn(u32 inst) {
	return (inst & (0b11111 << 5)) >> 5;
}

static Inst data_proc_imm(u32 binst) {
	Inst inst = UNKNOWN_INST;

	u32 op0 = (binst & (0b11 << 24)) >> 24;
	u32 op1 = (binst & (0b11 << 22)) >> 22;
	u32 top3 = (binst & (0b111 << 29)) >> 29;

	enum {
		Unknown,
		PCRelAddr,
		AddSub,
		LogicMove,
		BitfieldExtract
	} kind = Unknown;

	switch (op0) {
	case 0b00: kind = PCRelAddr;       break;
	case 0b01: kind = AddSub;          break;
	case 0b10: kind = LogicMove;       break;
	case 0b11: kind = BitfieldExtract; break;
	default:
		return UNKNOWN_INST;
	}

	switch (kind) {
	case PCRelAddr:
		if ((top3 & 0b100) == 0)
			inst.op = A64_ADR;
		else
			inst.op = A64_ADRP;

		inst.rd = regRd(binst);
		u64 immhi = binst & (0b111111111111111111 << 5);
		u64 immlo = binst & (0b11 << 29);
		inst.imm = (immhi >> 3) | (immlo >> 29);
		break;

	case AddSub:
		if (op1 == 0b10) {
			if (top3 == 0b100)
				inst.op = A64_ADDG;
			else if (top3 == 0b110)
				inst.op = A64_SUBG;
			else
				return UNKNOWN_INST;

			printf("ADDG, SUBG not supported"); // XXX
			return UNKNOWN_INST;
		}

		// normal ADDs and SUBs
		bool is_add = (top3 & 0b010) == 0;
		inst.op = (is_add) ? A64_ADD_IMM : A64_SUB_IMM;
		inst.w32 = (top3 & 0b100) == 0;
		inst.set_flags = (top3 & 0b001) > 0;

		u64 unshifted_imm = (binst & (0b111111111111 << 10)) >> 10;
		bool shift_by_12 = (binst & (1 << 22)) > 0;
		inst.imm = (shift_by_12) ? unshifted_imm << 12 : unshifted_imm;

		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
		break;

	case LogicMove:
		printf("Logic (imm) & move not supported"); break; // XXX

	case BitfieldExtract:
		printf("Bitfield & extract not supported"); break; // XXX
	}

	return inst;
}

// decode decodes n binary instructions from the input buffer and
// puts them in the output buffer, which must have space for n Insts.
int decode(u32 *in, int n, Inst *out) {
	int i;
	for (i = 0; i < n; i++) {
		u32 binst = in[i];
		switch (binst & (0b111 << 26)) {
		case 0b100 << 26:
			out[i] = data_proc_imm(binst);
			break;
		default:
			out[i] = UNKNOWN_INST;
		}
	}
	return i;
}

int main(int argc, char **argv) {
	int n = 1;
	u32 ibuf[] = {
		0b10010001000101010101010000100000, // ADD IMM X0, X1, #0b010101010101 (0x555)
	};
	Inst obuf[n];

	decode(ibuf, n, obuf);

	for (int i = 0; i < n; i++) {
		Inst inst = obuf[i];
		char regch = (inst.w32) ? 'W' : 'X';
		char flagsch = (inst.set_flags) ? 'S' : ' ';

		printf("%d%c %c%d, %c%d, #0x%llx\n", inst.op, flagsch,
			regch, inst.rd, regch, inst.rn, inst.imm);
	}

	return 0;
}
