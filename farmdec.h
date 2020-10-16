typedef enum AddrMode AddrMode;
typedef enum Cond Cond;
typedef enum ExtendType ExtendType;
typedef enum FPSize  FPSize;
typedef enum MemOrdering MemOrdering;
typedef enum Op          Op;
typedef enum Shift       Shift;
typedef enum Size Size;
typedef enum VectorArrangement VectorArrangement;
typedef struct Inst Inst;
typedef u8 Reg;

Cond get_cond(u8 flags);
AddrMode get_addrmode(u8 flags);
ExtendType get_mem_extend(u8 flags);

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
	A64_UNKNOWN, // unknown instruction (or Op field not set by accident)
	A64_ERROR,   // invalid instruction, Inst.error contains error string
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
	A64_BFM,
	A64_BFC,
	A64_BFI,
	A64_BFXIL,
	A64_UBFM,
	A64_LSL_IMM,
	A64_LSR_IMM,
	A64_UBFIZ,
	A64_UBFX,

	// Synthetic instruction comprising the SXTB, SXTH, SXTW, UXTB and UXTH aliases of SBFM and UBFM.
	// The kind of extension is stored in Inst.extend.type.
	A64_EXTEND,

	// Extract
	A64_EXTR,
	A64_ROR_IMM, // ROR Rd, Rs, #shift -- EXTR alias (Rm := Rs, Rn := Rs, predicate: Rm == Rn)

	/*** Branches, Exception Generating and System Instructions ***/

	A64_BCOND,

	// Exception generation
	//
	// With the exception of SVC, they are not interesting for lifting
	// userspace programs, but were included since they are trivial.
	A64_SVC, // system call
	A64_HVC,
	A64_SMC,
	A64_BRK,
	A64_HLT,
	A64_DCPS1,
	A64_DCPS2,
	A64_DCPS3,

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

	// Load-acquire/store-release register     -- AM_SIMPLE
	// Load/store register pair (post-indexed) -- AM_POST
	// Load/store register pair (offset)       -- AM_OFF_IMM
	// Load/store register pair (pre-indexed)  -- AM_PRE
	A64_LDP, // LDP, LDXP
	A64_STP, // STP, STXP
	A64_LDP_FP,
	A64_STP_FP,

	// Load/store register (unscaled immediate)
	A64_LDUR,
	A64_STUR,
	A64_LDUR_FP,
	A64_STUR_FP,

	// Load/store register (unprivileged): unsupported system instructions

	// Load register (literal)                      -- AM_LITERAL
	// Load-acquire/store-release register          -- AM_SIMPLE
	// Load-LOAcquire/Store-LORelease register      -- AM_SIMPLE
	// Load/store register (immediate post-indexed) -- AM_POST
	// Load/store register (immediate pre-indexed)  -- AM_PRE
	// Load/store register (register offset)        -- AM_OFF_REG, AM_OFF_EXT
	// Load/store register (unsigned immediate)     -- AM_OFF_IMM
	A64_LDR, // LDR, LDAR, LDLAR
	A64_STR, // STR, STLR, STLLR
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
	A64_LDSMAX,
	A64_LDSMIN,
	A64_LDUMAX,
	A64_LDUMIN,
	A64_SWP,
	A64_CAS,   // Compare and Swap (actually from Exclusive group)
	A64_CASP,  // Compare and Swap Pair of (double)words (actually from Exclusive group)
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



enum Shift {
	SH_LSL = 0b00,
	SH_LSR = 0b01,
	SH_ASR = 0b10,
	SH_ROR = 0b11, // only for RORV instruction; shifted add/sub does not support it
	SH_RESERVED = SH_ROR
};

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

enum flagmasks {
	W32 = 1 << 0,       // use the 32-bit W0...W31 facets?
	SET_FLAGS = 1 << 1, // modify the NZCV flags? (S mnemonic suffix)
};

// An Inst is a decoded instruction. Depending on the .op field, you can access
// the proper registers and immediates. Generally, if there is just one immediate,
// it is in the .imm field, or in the .offset field. If the instruction has more
// than one immediate, they are represented by one of the structs in the union.
struct Inst {
	Op op;

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
		}; // XXX according to the C standard, these three share their memory!
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
		char *error; // error string for op = A64_ERROR

		struct {
			u32 imm16;
			u32 lsl;   // left shift amount (0, 16, 32, 48)
		} mov_wide;
		struct {
			u32 lsb;
			u32 width;
		} bfm; // BFM aliases: BFXIL, SBFIZ, SBFX, UBFIZ, UBFX
		struct {
			u32 nzcv;
			u32 imm5;
		} ccmp;
		struct {
			s32 offset; // 14-bit jump offset
			u32 bit;    // b5:b40 field -- bit number to be tested
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

			Reg rs;    // status register for exclusive store (STXR, STLXR, STXP, STLXP)
		} ldst_order; // atomics and load/stores from Exclusive group
		struct {
			u32 nreg;   // consecutive vector registers to load/store
			u16 index;  // for single-struct variants: index of vector lane to load/store
			s16 offset; // offset to use if AM_POST and offset register Rm == ZERO_REG
		} simd_ldst; // LD1..4, ST1..4
	};
};
