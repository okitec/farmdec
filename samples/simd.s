simd_copy:
	dup d1, v2.d[1]
	dup s1, v2.s[3]
	dup v1.16b, v2.b[1]
	dup v1.8h, v2.h[3]
	dup v1.2s, v2.s[2]
	dup v1.2d, v2.d[0]
	dup v1.16b, w2
	dup v1.8h, w2
	dup v1.2s, w2
	dup v1.2d, x2
	ins v1.b[5], v2.b[3]
	ins v1.h[2], v2.h[1]
	ins v1.d[1], v2.d[0]
	ins v1.b[15], w2
	ins v1.h[3], w2
	ins v1.d[1], x2
	smov w1, v2.b[7]
	smov w1, v2.h[1]
	smov x1, v2.s[3]
	umov w1, v2.b[13]
	umov w1, v2.h[1]
	umov x1, v2.d[1]

threesame:
	shadd v1.2s, v2.2s, v3.2s
	uhadd v1.16b, v2.16b, v3.16b
	srhadd v1.2s, v2.2s, v3.2s
	urhadd v1.16b, v2.16b, v3.16b
	sqadd v1.2d, v2.2d, v3.2d
	sqadd s1, s2, s3
	uqadd v1.4h, v2.4h, v3.4h
	uqadd d1, d2, d3
	and v1.8b, v2.8b, v3.8b
	bic v1.8b, v2.8b, v3.8b
	orr v1.8b, v2.8b, v3.8b
	orn v1.8b, v2.8b, v3.8b
	eor v1.8b, v2.8b, v3.8b
	bsl v1.8b, v2.8b, v3.8b
	bit v1.8b, v2.8b, v3.8b
	bif v1.8b, v2.8b, v3.8b
	shsub v1.4s, v2.4s, v3.4s
	uhsub v1.8b, v2.8b, v3.8b
	sqsub v1.2d, v2.2d, v3.2d
	uqsub v1.8h, v2.8h, v3.8h
	cmgt v1.2s, v2.2s, v3.2s
	cmgt d1, d2, d3
	cmhi v1.2s, v2.2s, v3.2s
	cmhi d1, d2, d3
	cmge v1.2s, v2.2s, v3.2s
	cmge d1, d2, d3
	cmhs v1.2s, v2.2s, v3.2s
	cmhs d1, d2, d3
	sshl v1.2s, v2.2s, v3.2s
	sshl d1, d2, d3
	ushl v1.2s, v2.2s, v3.2s
	ushl d1, d2, d3
	sqshl v1.2s, v2.2s, v3.2s
	sqrshl v1.2s, v2.2s, v3.2s
	srshl v1.2s, v2.2s, v3.2s
	smax v1.2s, v2.2s, v3.2s
	umin v1.2s, v2.2s, v3.2s
	uabd v1.2s, v2.2s, v3.2s
	saba v1.2s, v2.2s, v3.2s
	add v1.2s, v2.2s, v3.2s
	add d1, d2, d3
	sub v1.2s, v2.2s, v3.2s
	cmtst v1.2s, v2.2s, v3.2s
	cmeq v1.2s, v2.2s, v3.2s
	mla v1.2s, v2.2s, v3.2s
	mls v1.2s, v2.2s, v3.2s
	mul v1.2s, v2.2s, v3.2s
	pmul v1.16b, v2.16b, v3.16b
	smaxp v1.2s, v2.2s, v3.2s
	uminp v1.2s, v2.2s, v3.2s
	sqdmulh v1.2s, v2.2s, v3.2s
	sqrdmulh v1.2s, v2.2s, v3.2s
	addp v1.2s, v2.2s, v3.2s
	fmaxnm v1.2s, v2.2s, v3.2s
	fminnm v1.2s, v2.2s, v3.2s
	fmaxnmp v1.2s, v2.2s, v3.2s
	fmla v1.2s, v2.2s, v3.2s
	fmls v1.2s, v2.2s, v3.2s
	fmlal v1.2s, v2.2h, v3.2h
	fmlsl v1.2s, v2.2h, v3.2h
	fadd v1.2s, v2.2s, v3.2s
	fsub v1.2s, v2.2s, v3.2s
	fabd v1.2s, v2.2s, v3.2s
	fmulx v1.2s, v2.2s, v3.2s
	fmulx d1, d2, d3
	fmul v1.2s, v2.2s, v3.2s
	fcmeq v1.2s, v2.2s, v3.2s
	fcmge v1.2s, v2.2s, v3.2s
	fcmgt v1.2s, v2.2s, v3.2s
	facge v1.2s, v2.2s, v3.2s
	facgt v1.2s, v2.2s, v3.2s
	fmax v1.2s, v2.2s, v3.2s
	fmin v1.2s, v2.2s, v3.2s
	frecps v1.2s, v2.2s, v3.2s
	frsqrts v1.2s, v2.2s, v3.2s
	fdiv v1.2s, v2.2s, v3.2s

extract:
	ext v1.8b, v2.8b, v3.8b, #7
	ext v1.16b, v2.16b, v3.16b, #15

table:
	tbl v1.8b, {v2.16b}, v3.8b
	tbl v1.16b, {v2.16b}, v3.16b
	tbl v1.8b, {v2.16b, v3.16b}, v4.8b
	tbl v1.8b, {v2.16b, v3.16b, v4.16b}, v5.8b
	tbl v1.8b, {v2.16b, v3.16b, v4.16b, v5.16b}, v6.8b
	tbx v1.8b, {v2.16b}, v3.8b
	tbx v1.16b, {v2.16b}, v3.16b
	tbx v1.8b, {v2.16b, v3.16b}, v4.8b
	tbx v1.8b, {v2.16b, v3.16b, v4.16b}, v5.8b
	tbx v1.8b, {v2.16b, v3.16b, v4.16b, v5.16b}, v6.8b

permute:
	uzp1 v1.2d, v2.2d, v3.2d
	uzp2 v1.2s, v2.2s, v3.2s
	trn1 v1.8h, v2.8h, v3.8h
	trn2 v1.4h, v2.4h, v3.4h
	zip1 v1.16b, v2.16b, v3.16b
	zip2 v1.8b, v2.8b, v3.8b

threesameextra:
	sqrdmlah v1.4s, v2.4s, v3.4s
	sqrdmlah h1, h2, h2
	sqrdmlsh v1.8h, v2.8h, v3.8h
	sqrdmlsh s1, s2, s3
	sdot v1.2s, v2.8b, v3.8b
	udot v1.4s, v2.16b, v3.16b
	fcmla v1.2d, v2.2d, v3.2d, #0
	fcmla v1.2d, v2.2d, v3.2d, #90
	fcmla v1.2d, v2.2d, v3.2d, #180
	fcmla v1.2d, v2.2d, v3.2d, #270
	fcadd v1.4s, v2.4s, v3.4s, #90
	fcadd v1.4s, v2.4s, v3.4s, #270
