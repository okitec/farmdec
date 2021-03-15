#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define FARMDEC_INTERNAL
#include "farmdec.h"

enum {
	MAX_OBUF_SIZE = 1 * 1024 * 1024 * 1024 // 1 GB
};

// char *mnemonics[], generated using mnemonics.awk.
#include "mnemonics.h"

void printinst(Inst inst);

// main: read instructions as list of hexadecimal numbers, decode, print results.
int main(int argc, char **argv) {
	// Single instruction as argument?
	if (argc == 2) {
		u32 binst;
		Inst inst;
		sscanf(argv[1], "%x", &binst);
		fad_decode(&binst, 1, &inst);
		printinst(inst);
		return 0;
	}

	// ----- Benchmark/Bulk mode -----

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

	fad_decode(ibuf, ninst, obuf);

	for (uint i = 0; i < ninst; i++) {
		printf("%06x ", 4*i); // addr in ibuf
		printinst(obuf[i]);
	}

	return 0;
}

void printinst(Inst inst) {
	char regch = (inst.flags & W32) ? 'W' : 'X';
	char flagsch = (inst.flags & SET_FLAGS) ? 'S' : ' ';

	switch (inst.op) {
	case A64_UNKNOWN: {
		u32 binst = inst.imm;

		// Print unknown instruction as little-endian hex bytes. This format
		// can be used by the Online Disassembler (ODA).
		printf("??? %02x %02x %02x %02x\n", (binst>>0) & 0xff, (binst>>8) & 0xff,
			(binst>>16) & 0xff, (binst>>24) & 0xff);
		return;
	}
	case A64_ERROR:   printf("error \"%s\"\n", inst.error); return;
	case A64_UDF:     printf("udf\n");                      return;
	default:
		break; // normal instruction
	}

	// We do not disambiguate here -- all instructions are printed
	// the same; for example, instructions with two immediates have
	// the imm field printed too.
	printf("%-12s %c %c%d, %c%d, %c%d, imm=%lu, offset=%+ld, fimm=%f, imm2=(%u,%u), flags=0b %o%o%o %o%o%o %o%o\n",
		mnemonics[inst.op], flagsch,
		regch, inst.rd, regch, inst.rn, regch, inst.rm,
		inst.imm, inst.offset, inst.fimm, inst.bfm.lsb, inst.bfm.width,
		(inst.flags>>7) & 1, (inst.flags>>6) & 1, (inst.flags>>5) & 1, (inst.flags>>4) & 1,
		(inst.flags>>3) & 1, (inst.flags>>2) & 1, (inst.flags>>1) & 1, (inst.flags>>0) & 1);
}
