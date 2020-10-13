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

typedef enum AddrMode AddrMode;
typedef enum Cond Cond;
typedef enum ExtendType ExtendType;
typedef enum FPSize FPSize;
typedef enum Op Op;
typedef enum Shift Shift;
typedef enum Size Size;
typedef enum VectorArrangement VectorArrangement;
typedef struct Inst Inst;
typedef u8 Reg;

// Opcodes ordered and grouped according to the Top-level Encodings
// of the A64 Instruction Set Architecture (ARMv8-A profile) document,
// pages 1406-1473.
//
// Immediate and register variants generally have different opcodes
// (e.g. A64_ADD_IMM, A64_ADD_SHIFTED, A64_ADD_EXT), but the marker
// only appears where disambiguation is needed (e.g. ADR is not called
// ADR_IMM since there is no register variant). Aliases have an opcode
// of their own.
//
// Where possible, variants of instructions with regular structure
// are encoded as one instruction. For example, conditional branches
// like B.EQ, B.PL and so on are encoded as A64_BCOND with the
// condition encoded in the Inst.flags field. The various addressing
// modes of loads and stores are encoded similarly. See the Inst
// structure for more detail.
enum Op {
	A64_UNKNOWN, // an invalid (Inst.error != NULL) or (to us) unknown instruction
	A64_UDF,     // throws undefined exception

	/*** Data Processing -- Immediate ***/

	// PC-rel. addressing
	A64_ADR,     // ADR Xd, label  -- Xd ← PC + label
	A64_ADRP,    // ADRP Xd, label -- Xd ← PC + (label * 4K)

	// Add/subtract (immediate, with tags) -- OMITTED

	// Add/subtract (immediate)
	A64_ADD_IMM,
	A64_MOV_SP, // MOV from/to SP -- ADD (imm) alias (predicate: shift == 0 && imm12 == 0 && (Rd == SP || Rn == SP))
	A64_SUB_IMM,

	// Logical (immediate)
	A64_AND_IMM,
	A64_ORR_IMM,
	A64_EOR_IMM,
	A64_TST_IMM, // TST Rn -- ANDS alias (Rd := RZR, predicate: Rd == ZR && set_flags)

	// Move wide (immediate)
	A64_MOVN,
	A64_MOVZ,
	A64_MOVK,

	// Synthetic instruction comprising MOV (bitmask immediate), MOV (inverted wide immediate)
	// and MOV (wide immediate), which are themselves aliases of ORR_IMM, MOVN and MOVZ respectively.
	// For lifting, we do not care how the immediate is encoded, only that it is an immediate move.
	A64_MOV_IMM,

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
	A64_ROR_IMM, // ROR Rd, Rs, #shift -- EXTR alias (Rm := Rs, Rn := Rs, predicate: Rm == Rn)

	/*** Branches, Exception Generating and System Instructions ***/

	A64_BCOND,

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
	A64_TBNZ,

	/*** Data Processing -- Register ***/

	// Data-processing (2 source)
	A64_UDIV,
	A64_SDIV,
	A64_LSLV,
	A64_LSRV,
	A64_ASRV,
	A64_RORV,
	A64_CRC32B,
	A64_CRC32H,
	A64_CRC32W,
	A64_CRC32X,
	A64_CRC32CB,
	A64_CRC32CH,
	A64_CRC32CW,
	A64_CRC32CX,
	A64_SUBP,

	// Data-processing (1 source)
	A64_RBIT,
	A64_REV16,
	A64_REV,
	A64_REV32,
	A64_CLZ,
	A64_CLS,

	// Logical (shifted register)
	A64_AND_SHIFTED,
	A64_TST_SHIFTED, // ANDS alias (Rd := ZR, predicate: Rd == ZR)
	A64_BIC,
	A64_ORR_SHIFTED,
	A64_MOV_REG,     // ORR alias (predicate: shift == 0 && imm6 == 0 && Rn == ZR)
	A64_ORN,
	A64_MVN,         // ORN alias (Rn := ZR, predicate: Rn == ZR)
	A64_EOR_SHIFTED,
	A64_EON,

	// Add/subtract (shifted register)
	A64_ADD_SHIFTED,
	A64_CMN_SHIFTED, // ADDS alias (Rd := ZR, predicate: Rd == ZR && set_flags)
	A64_SUB_SHIFTED,
	A64_NEG,         // SUB alias (Rn := ZR, predicate: Rn == ZR)
	A64_CMP_SHIFTED, // SUBS alias (Rd := ZR, predicate: Rd == ZR && set_flags)

	// Add/subtract (extended register)
	// Register 31 is interpreted as the stack pointer (SP/WSP).
	A64_ADD_EXT,
	A64_CMN_EXT, // ADDS alias (Rd := ZR, predicate: Rd == ZR && set_flags)
	A64_SUB_EXT,
	A64_CMP_EXT, // SUBS alias (Rd := ZR, predicate: Rd == ZR && set_flags)

	// Add/subtract (with carry)
	A64_ADC,
	A64_SBC,
	A64_NGC, // SBC alias (Rd := ZR, predicate: Rd == RR)

	// Rotate right into flags
	A64_RMIF,

	// Evaluate into flags
	A64_SETF8,
	A64_SETF16,

	// Conditional compare (register)
	A64_CCMN_REG,
	A64_CCMP_REG,

	// Conditional compare (immediate)
	A64_CCMN_IMM,
	A64_CCMP_IMM,

	// Conditional select
	A64_CSEL,
	A64_CSINC,
	A64_CINC,  // CSINC alias (cond := invert(cond), predicate: Rm == Rn != ZR)
	A64_CSET,  // CSINC alias (cond := invert(cond), predicate: Rm == Rn == ZR)
	A64_CSINV,
	A64_CSETM, // CSINV alias (cond := invert(cond), predicate: Rm == Rn == ZR)
	A64_CSNEG,
	A64_CNEG,  // CSNEG alias (cond := invert(cond), predicate: Rm == Rn)

	// Data-processing (3 source)
	A64_MADD,
	A64_MUL,    // MADD alias (Ra omitted, predicate: Ra == ZR)
	A64_MSUB,
	A64_MNEG,   // MSUB alias (^---- see above)
	A64_SMADDL,
	A64_SMULL,  // SMADDL alias  (^---- see above)
	A64_SMSUBL,
	A64_SMNEGL, // SMSUBL alias (^---- see above)
	A64_SMULH,
	A64_UMADDL,
	A64_UMULL,  // UMADDL alias (^---- see above)
	A64_UMSUBL,
	A64_UMNEGL, // UMSUBL alias (^---- see above)
	A64_UMULH,

	/*** Loads and Stores ***/

	// There are not that many opcodes because access size, sign-extension
	// and addressing mode (post-indexed, register offset, immediate) are
	// encoded in the Inst, to leverage the regular structure and cut down
	// on opcodes (and by extension, duplicative switch-cases for the user
	// of this decoder).

	// Advanced SIMD load/store multiple structures
	// Advanced SIMD load/store multiple structures (post-indexed)
	A64_LD1_MULT,
	A64_ST1_MULT,
	A64_LD2_MULT,
	A64_ST2_MULT,
	A64_LD3_MULT,
	A64_ST3_MULT,
	A64_LD4_MULT,
	A64_ST4_MULT,

	// Advanced SIMD load/store single structure
	// Advanced SIMD load/store single structure (post-indexed)
	A64_LD1_SINGLE,
	A64_ST1_SINGLE,
	A64_LD2_SINGLE,
	A64_ST2_SINGLE,
	A64_LD3_SINGLE,
	A64_ST3_SINGLE,
	A64_LD4_SINGLE,
	A64_ST4_SINGLE,
	A64_LD1R,
	A64_LD2R,
	A64_LD3R,
	A64_LD4R,

	// Load/store exclusive
	A64_LDXR,  // includes Load-acquire variants
	A64_STXR,  // includes Store-acquire variants (STLXR)
	A64_LDXP,  // ------
	A64_STXP,  // ------
	A64_LDAPR, // Load-AcquirePC Register (actually in Atomic group)

	// Load/store no-allocate pair (offset)
	A64_LDNP,
	A64_STNP,
	A64_LDNP_FP,
	A64_STNP_FP,

	// Load/store register pair (post-indexed)
	// Load/store register pair (offset)
	// Load/store register pair (pre-indexed)
	A64_LDP,
	A64_STP,
	A64_LDP_FP,
	A64_STP_FP,

	// Load/store register (unscaled immediate)
	A64_LDUR,
	A64_STUR,
	A64_LDUR_FP,
	A64_STUR_FP,

	// Load/store register (unprivileged): unsupported system instructions

	// Load register (literal)
	// Load-acquire/store-release register
	// Load-LOAcquire/Store-LORelease register
	// Load/store register (immediate post-indexed)
	// Load/store register (immediate pre-indexed)
	// Load/store register (register offset)
	// Load/store register (unsigned immediate)
	A64_LDR,
	A64_STR,
	A64_LDR_FP,
	A64_STR_FP,

	// Atomic memory operations
	//
	// Whether the instruction has load-acquire (e.g. LDADDA*), load-acquire/
	// store-release (e.g. LDADDAL*) or store-release (e.g. STADDL) semantics
	// is stored in ldst_order.load and .store.
	//
	// There are no ST* aliases; the only difference to the LD* instructions
	// is that the original value of the memory cell is discarded by writing
	// to the zero register.
	A64_LDADD,
	A64_LDCLR,
	A64_LDEOR,
	A64_LDSET,
	A64_LDSMAX, // Signed/Unsigned (SMAX, UMAX) → flags, not opcode ? XXX
	A64_LDSMIN,
	A64_LDUMAX,
	A64_LDUMIN,
	A64_SWP,
	A64_CAS,   // Compare and Swap (actually from Exclusive group)
	A64_CASP,  // Compare and Swap Pair of (double)words (actually from Exclusive group)
};

enum Shift {
	SH_LSL = 0b00,
	SH_LSR = 0b01,
	SH_ASR = 0b10,
	SH_ROR = 0b11, // only for RORV instruction; shifted add/sub does not support it
	SH_RESERVED = SH_ROR
};

// Size, encoded in two bits.
enum Size {
	SZ_B = 0b00, // Byte     -  8 bit
	SZ_H = 0b01, // Halfword - 16 bit
	SZ_W = 0b10, // Word     - 32 bit
	SZ_X = 0b11, // Extended - 64 bit
};

// Floating-point size, encoded in three bits. Mostly synonymous to Size, but
// with the 128-bit quadruple precision.
enum FPSize {
	FSZ_B = SZ_B, // Byte   -   8 bits
	FSZ_H = SZ_H, // Half   -  16 bits
	FSZ_S = SZ_W, // Single -  32 bits
	FSZ_D = SZ_X, // Double -  64 bits

	// "Virtual" encoding, never used in the actual instructions.
	// There, Quad precision is encoded in various incoherent ways.
	FSZ_Q = 0b111 // Quad   - 128 bits
};

// The three-bit Vector Arrangement specifier determines the structure of the
// vectors used by a SIMD instruction, where it is encoded in size(2):Q(1).
//
// The vector registers V0...V31 are 128 bit long, but some arrangements use
// only the bottom 64 bits.
enum VectorArrangement {
	VA_8B  = (FSZ_B << 1) | 0, //  64 bit
	VA_16B = (FSZ_B << 1) | 1, // 128 bit
	VA_4H  = (FSZ_H << 1) | 0, //  64 bit
	VA_8H  = (FSZ_H << 1) | 1, // 128 bit
	VA_2S  = (FSZ_S << 1) | 0, //  64 bit
	VA_4S  = (FSZ_S << 1) | 1, // 128 bit
	VA_1D  = (FSZ_D << 1) | 0, //  64 bit
	VA_2D  = (FSZ_D << 1) | 1, // 128 bit
};

// ExtendType: signed(1):size(2)
enum ExtendType {
	UXTB = (0 << 2) | SZ_B,
	UXTH = (0 << 2) | SZ_H,
	UXTW = (0 << 2) | SZ_W,
	UXTX = (0 << 2) | SZ_X,
	SXTB = (1 << 2) | SZ_B,
	SXTH = (1 << 2) | SZ_H,
	SXTW = (1 << 2) | SZ_W,
	SXTX = (1 << 2) | SZ_X,
};

// XXX keep it at 16 bytes if possible
struct Inst {
	Op  op;

	// Overloaded flags bitfield. The two lowest bits W32 and SET_FLAGS are never overloaded.
	//
	//              7   6   5   4   3   2   1   0
	// Default:     - | - | - | - | - | - | S | W32
	// Conditional: -----cond-----| - | - | S | W32  (see enum Cond)
	// Load/Store:  ---mode---|----ext----|N/A| W32  (see enum AddrMode, enum ExtendType)
	//     The load/store ext field stores the access size and and whether to do sign
	//     extension or zero extension. Hence it is identical to an ExtendType.
	// L/S Float:   ---mode---|----size---|N/A| N/A  (see enum FPSize)
	//     The SIMD+FP variants of the usual LDR, STR, ...
	// LDx/STx:     ---mode---|----vec----|N/A| N/A  (see enum VectorArrangement)
	//     LD1..4, ST1..4, etc., need the vector arrangement, like many SIMD operations.
	u8 flags;

	union {
		struct {
			Reg rd;   // destination register - Rd
			Reg rn;   // first (or only) operand, read-only - Rn, Rt (CBZ)
			Reg rm;   // second operand, read-only - Rm
		};
		struct {
			Reg rt;  // destination of load, source of store
			Reg rn;  // base addressing register (Xn)
			union {
				Reg rt2; // second destination/source register for LDP, STP and variants (e.g. LDXP)
				Reg rm;  // index register for AM_OFF_REG, AM_OFF_EXT
				Reg rs;  // source register for atomic operations
			};
		} ldst;
	};
	union {
		u64 imm;     // single immediate
		s64 offset;  // branches, ADR, ADRP: PC-relative byte offset
		Reg ra;      // third operand for 3-source data proc instrs (MADD, etc.)
		char *error; // error string for op = A64_UNKNOWN, may be NULL

		struct {
			u32 imm16;
			u32 lsl;   // left shift amount (0, 16, 32, 48)
		} mov_wide;
		struct {
			u32 immr;
			u32 imms;
		}; // XXX just used by bitfield operations, right?
		struct {
			u32 nzcv;
			u32 imm5;
		} ccmp;
		struct {
			s32 offset; // 14-bit jump offset
			u32 b5b40;  // b5:b40 field - bit number to be tested
		} tbz;
		struct {
			u32 type;   // enum Shift (not used because sizeof(enum) is impl-defined)
			u32 amount;
		} shift;
		struct {
			u32 mask;
			u32 ror;  // rotate right amount - 0..63
		} rmif;
		struct {
			u32 type; // enum ExtendType
			u32 lsl;  // left shift amount
		} extend;
		struct {
			// Atomics can have different ordering for the load and the store, so
			// we need to have two variables.
			u16 load;  // enum MemOrdering for Load -- None, Acquire, LOAcquire, AcquirePC
			u16 store; // --------------- for Store -- None, Release, LORelease

			Reg rs;    // status register for exclusive store pair (STXP, STLXP)
		} ldst_order; // atomics and load/stores from Exclusive group
		struct {
			u32 nreg;   // consecutive vector registers to load/store
			u16 index;  // for single-struct variants: index of vector lane to load/store
			s16 offset; // offset to use if AM_POST and offset register Rm == ZERO_REG
		} simd_ldst; // LD1..4, ST1..4
	};
};

static Inst UNKNOWN_INST = {
	op: A64_UNKNOWN,
	// all other fields: zero
};

static Inst errinst(char *err) {
	Inst inst = UNKNOWN_INST;
	inst.error = err;
	return inst;
}

enum flagmasks {
	W32 = 1 << 0,       // use the 32-bit W0...W31 facets?
	SET_FLAGS = 1 << 1, // modify the NZCV flags? (S mnemonic suffix)
};

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
	cond &= 0xF;
	flags &= 0x0F;
	return (cond << 4) | flags;
}

static u8 invert_cond(u8 flags) {
	Cond cond = get_cond(flags);
	return set_cond(flags, cond ^ 0b001); // invert LSB
}

// Addressing modes, stored in the top three bits of the flags field
// (where the condition is stored for conditional instructions). See
// page 187, section C1.3.3 of the 2020 ARM ARM for ARMv8.
//
// The base register is stored in the Inst.ldst.rn field.
//
// The LSL amount for the REG and EXT depends on the access size
// (#4 for 128 bits (SIMD), #3 for 64 bits, #2 for 32 bits, #1 for
// 16 bits, #0 for 8 bits) and is used for array indexing:
//
//     u64 a[128];
//     u64 x0 = a[i]; → ldr x0, [a, i, LSL #3]
//
enum AddrMode {
	AM_SIMPLE,  // [base] -- used by atomics, exclusive, ordered load/stores → check Inst.ldst_order
	AM_OFF_IMM, // [base, #imm]
	AM_OFF_REG, // [base, Xm, {LSL #imm}] (#imm either #log2(size) or #0)
	AM_OFF_EXT, // [base, Wm, {S|U}XTW {#imm}] (#imm either #log2(size) or #0)
	AM_PRE,     // [base, #imm]!
	AM_POST,    // [base],#imm  (for LDx, STx also register: [base],Xm)
	AM_LITERAL  // label
};

// Addressing mode, for Loads and Stores.
static AddrMode get_addrmode(u8 flags) {
	return (AddrMode)((flags >> 5) & 0b111);
}

static u8 set_addrmode(u8 flags, AddrMode mode) {
	return ((mode&0b111) << 5) | (flags&0b11111);
}

// How much memory to load/store (access size) and whether to sign- or
// zero-extend the value.
static ExtendType get_mem_extend(u8 flags) {
	return (ExtendType)((flags >> 2) & 0b111);
}

static u8 set_mem_extend(u8 flags, ExtendType memext) {
	return ((memext&0b111) << 2) | (flags&0b11100011);
}

// Memory ordering semantics for Atomic instructions and the Load/Stores in the
// Exclusive group.
enum MemOrdering {
	MO_NONE,
	MO_ACQUIRE,    // Load-Acquire -- sequentially consistent Acquire
	MO_LO_ACQUIRE, // Load-LOAcquire -- Acquire in Limited Ordering Region (LORegion)
	MO_ACQUIRE_PC, // Load-AcquirePC -- weaker processor consistent (PC) Acquire
	MO_RELEASE,    // Store-Release
	MO_LO_RELEASE, // Store-LORelease -- Release in LORegion
};

// Register 31's interpretation is up to the instruction. Many interpret it as the
// zero register ZR/WZR. Reading to it yields a zero, writing discards the result.
// Other instructions interpret it as the stack pointer SP.
//
// We split up this overloaded register: when we encounter R31 and interpret it as
// the stack pointer, we assign a different number. This way, the user does not
// need to know which instructions use the SP and which use the ZR.
enum special_registers {
	ZERO_REG = 31,
	STACK_POINTER = 100 // arbitrary
};


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
			default:
				// Shut up impossible unhandled enum value warning.
				return errinst("data_proc_imm/Move: cannot happen");
			}
		}

		break;
	}
	case Bitfield: {
		switch (top3 & 0b011) { // base instructions at this point
		case 0b00: inst.op = A64_SBFM; break;
		case 0b01: inst.op = A64_BFM;  break;
		case 0b10: inst.op = A64_UBFM; break;
		default:
			return errinst("data_proc_imm/Bitfield: neither SBFM, BFM or UBFM");
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
				inst.immr = 0xDEAD;  // XXX immr = -lsb MOD 32/64 → lsb = ?
				inst.imms = imms + 1; // imms = width - 1 → width = imms + 1
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
			inst.immr = immr;            // immr = lsb
			inst.imms = imms - immr + 1; // imms = lsb + width - 1 → width = imms - lsb + 1;

			break;
			// There is no SBFM general case not handled by an alias
			// XXX is that true?
		case A64_BFM:
			// XXX BFM aliases not implemented
			break;

		case A64_UBFM:
			// XXX UBFM aliases not implemented
			// XXX how to re-use parts of SBFM?
			break;
		default:
			return errinst("data_proc_imm/Bitfield: cannot happen");
		}

	regs:
		inst.rd = regRd(binst);
		inst.rn = regRn(binst);
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

	case System:
		return errinst("System instructions not supported"); // XXX _some_ should be decoded. barriers, CFINV

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
		inst.tbz.b5b40 = (b5 >> (31-5)) | b40; // b5:b40
		inst.rn = binst & 0b11111;             // Rt; not modified → Inst.rn, not .rd
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
			out[i] = errinst("Loads & stores not supported");
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

		if (inst.op == A64_UNKNOWN) {
			if (inst.error != NULL) {
				printf("0x%04x: error \"%s\"\n", 4*i, inst.error);
			} else {
				printf("0x%04x: ???\n", 4*i);
			}
			continue;
		}

		if (inst.op == A64_UDF) {
			printf("%04x udf\n", 4*i);
			continue;
		}

		// We do not disambiguate here -- all instructions are printed
		// the same; for example, instructions with two immediates have
		// the imm field printed too.
		printf("%04x %-12s%c %c%d, %c%d, %c%d, imm=%lu, imm2=(%u,%u), offset=%04ld, w32=%o, set_flags=%o, memext=%o, addrmode=%o, cond=%x\n",
			4*i, mnemonics[inst.op], flagsch,
			regch, inst.rd, regch, inst.rn, regch, inst.rm,
			inst.imm, inst.immr, inst.imms, inst.offset, inst.flags&W32, inst.flags&SET_FLAGS,
			get_mem_extend(inst.flags), get_addrmode(inst.flags), get_cond(inst.flags));
	}

	return 0;
}
