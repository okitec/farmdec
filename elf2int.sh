#!/bin/sh
#
# elf2int.sh <ELF file> -- convert ELF .text section of little-endian A64 instructions
# to list of unsigned integers that can be read using scanf("%x"). 

convert='
$1 ~ /0x.+/ {
	print from_little_endian($2)
	print from_little_endian($3)
	print from_little_endian($4)
	print from_little_endian($5)
}

# ff8300d1 → d18300ff
function from_little_endian(s,            a) {
	split(s, a, "")                                # split into single characters
	return a[7] a[8] a[5] a[6] a[3] a[4] a[1] a[2] # awk is 1-based
}
'

readelf -x .text $1 | awk "$convert"
