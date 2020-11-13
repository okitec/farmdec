# farmdec -- Fast ARM A64 Decoder

This decoder aims to be the simplest possible, without dependencies or
a complex API; its purpose is to be the A64 equivalent of the x86-64
decoder fadec, that is, to be used by the Rellume lifter.

### Usage

The `fad_decode` function takes an input and an output buffer. The input
is a buffer of binary A64 instructions as u32 (== uint32_t) and is hence
endianness-agnostic. If the caller receives the instructions as a byte
buffer, it must convert the little-endian bytes to the native endianness
in the u32 buffer.

The decoded instructions are of type Inst, which fits in 16 bytes and is
heavily overloaded; you need to look at the opcode first. The opcode should
give you enough clue which fields of the Inst to access (an A64_B branch
needs the `offset` field, an A64_ADD_IMM needs the `imm` immediate, an
A64_AND_REG needs the three registers). Refer to the extensive comments.

### Structure

The actual code is large, but extremely regular; it is essentially a shallow
hierarchical lookup table made of code, comparable to a recursive descent
parser. It follows the structure of the section "Top-level encodings for A64"
at the end of the "ARM A64 Instruction Set Architecture" document.

As an example, we shall decode this:

	ADD X0, X4, #32

1. We determine the instruction category ("Data Processing -- Immediate")
   using a switch on the bit patterns given by the document, and call
   the `data_proc_imm` function.
2. The function determines the kind of the instruction ("Add/subtract
   (immediate)") by another switch on bit patterns.
3. Given the kind, the exact opcode and all fields are determined (here:
   A64_ADD_IMM, with imm=32, rn=4, rd=0 and no flags set). If there are
   alias instructions, we select the preferred one based on the conditions
   given in the document.
4. The `data_proc_imm` returns the Inst to `decode`, which puts it in the
   output buffer.

The bit patterns in the document often have wildcards denoted with an `x`,
where the bit is variable and the value is not important at that stage of
decoding. For example, the pattern `100x` describes Immediate Data Processing
instructions while `x1x0` describes Loads and Stores. Instead of checking
single bits in if-else ladders, we construct switches with as-long-as-possible
bit patterns for which we can enumerate the individual values without much
trouble (one or two variable bits). Checking for `x1x0` works like this:

	case 0b0100:
	case 0b0110:
	case 0b1100:
	case 0b1110:
		do_something();
		break;

#### Rationale

The brute-force approach is simple to read and maintain even for people
unfamiliar with C -- there are no clever parts: not in the code itself,
not in the build process. Furthermore, the specification document can
serve as the documentation for this. The performance is also adequate,
with 10 ns per instruction on average.

The LLVM project has a lot of code for ARM A64, including an assembler
and a disassembler, but it is a lot of code with a much more complex
interface.

### Omissions

Some instructions are deliberately not handled due to a focus on
userspace programs.

 - Instructions working with memory tags (e.g. ADDG).
 - ERET
 - Instructions with Pointer Authentication (PAC) (e.g. AUTDA, PACIA).
 - Unprivileged load/store, which is used by OS kernels and supervisors.

Support for SIMD and Floating Point Instructions is missing but
will be added.
