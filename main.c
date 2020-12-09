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

	fad_decode(ibuf, ninst, obuf);

	for (uint i = 0; i < ninst; i++) {
		Inst inst = obuf[i];
		char regch = (inst.flags & W32) ? 'W' : 'X';
		char flagsch = (inst.flags & SET_FLAGS) ? 'S' : ' ';


		switch (inst.op) {
		case A64_UNKNOWN: printf("%04x ???\n", 4*i);                      continue;
		case A64_ERROR:   printf("%04x error \"%s\"\n", 4*i, inst.error); continue;
		case A64_UDF:     printf("%04x udf\n", 4*i);                      continue;
		default:
			break; // normal instruction
		}

		// We do not disambiguate here -- all instructions are printed
		// the same; for example, instructions with two immediates have
		// the imm field printed too.
		printf("%04x %-12s %c %c%d, %c%d, %c%d, imm=%lu, fimm=%f, imm2=(%u,%u), offset=%+ld, w32=%o, set_flags=%o, memext=%o, fpprec=%o, addrmode=%o, cond=%x\n",
			4*i, mnemonics[inst.op], flagsch,
			regch, inst.rd, regch, inst.rn, regch, inst.rm,
			inst.imm, inst.fimm, inst.bfm.lsb, inst.bfm.width, inst.offset, inst.flags&W32, inst.flags&SET_FLAGS,
			fad_get_mem_extend(inst.flags), fad_get_prec(inst.flags), fad_get_addrmode(inst.flags), fad_get_cond(inst.flags));
	}

	return 0;
}
