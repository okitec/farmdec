#!/usr/bin/awk -f
#
# mnemonics.awk <farmdec.h> -- print mnemonics string table for A64 Ops.

BEGIN {
	print "char *mnemonics[] = {"
}

$1 ~ /A64_.+/ {
	sub(/A64_/, "", $1)
	sub(/,/, "", $1)
	printf "\t\"%s\",\n", tolower($1)
}

END {
	print "};"
}
