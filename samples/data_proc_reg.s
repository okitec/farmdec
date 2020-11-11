two_source:
	udiv xzr, x1, x2
	sdiv xzr, x1, x2
	lslv wzr, w1, w2
	lsrv xzr, x1, x2
	asrv x10, x1, x2
	rorv xzr, x1, x2
	crc32b w3, w1, w2
	crc32h w3, w1, w2
	crc32w w3, w1, w2
	crc32x w3, w1, x2
	crc32cb wzr, w1, w2
	crc32ch w3, w1, w2
	crc32cw w3, w1, w2
	crc32cx wzr, w1, x2
	subp x10, x13, sp
	subps x10, sp, sp

one_source:
	rbit w1, w2
	rev16 x1, x2
	rev x1, x2
	clz x1, x2
	cls x1, x2
	rev32 x1, x3

logical:
	and w1, w2, w3
	bic x1, x2, x3, lsl #17
	orr x1, x2, x3, lsr #2
	orn x1, x2, x3, asr #63
	eor x1, x2, x3, ror #2
	eon wzr, w2, w3
	ands x1, x2, x3
	bics x1, x2, x3

add_sub_shifted:
	add xzr, x2, x3
	adds x1, x2, x3, lsr #3
	cmn x2, x3, lsr #4
	sub w1, w2, w3
	subs x1, x2, xzr
	cmp x2, x3

add_sub_extended:
	add sp, sp, x3
	adds x1, x2, w3, sxth
	cmn x2, w3, uxtb
	sub x1, sp, x3, sxtx
	subs x1, x2, w3, uxtb
	cmp x1, w2, sxtw

add_sub_carry:
	adc w1, w2, w3
	adcs x1, x2, x3
	sbc w1, w2, w3
	sbcs x1, x2, x3
	ngc w1, w2
	ngcs x1, x2

flags:
	rmif x3, #14, #3
	setf8 w4
	setf16 w4

cond_comp:
	ccmn x2, x3, #2, ne
	ccmp x2, x3, #12, ge

cond_comp_imm:
	ccmn x1, #12, #2, mi
	ccmp w1, #22, #4, pl

cond_sel:
	csel x1, x2, x3, lo
	csinc w1, w2, w3, vs
	cinc x1, x2, hi
	cset w1, cc
	csinv x1, x2, x3, vc
	cinv w1, w2, ls
	csetm x1, lt
	csneg w1, w2, w3, ge
	cneg x1, x2, hs

three_source:
	madd x1, x2, x3, x4
	mul x1, x2, x3
	msub w1, w2, w3, w4
	mneg x1, x2, x3
	smaddl x1, w2, w3, x4
	smull x1, w2, w3
	smsubl x1, w2, w3, x4
	smnegl x1, w2, w3
	smulh x1, x2, x3
	umaddl x1, w2, w3, x4
	umsubl x1, w2, w3, x4
	umulh x1, x2, x3
	umull x1, w2, w3
	umnegl x1, w2, w3

