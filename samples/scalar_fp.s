fmov:
	fmov h1, w2
	fmov s1, w2
	fmov d1, x2
	fmov w2, h1
	fmov w2, s1
	fmov x2, d1
	fmov v1.d[1], x2
	fmov x2, v1.d[1]

fp_fixed_point_conv:
	fcvtzs x2, d1, #15
	fcvtzu w2, h1, #3
	scvtf s1, x2, #17
	ucvtf d1, w2, #5

int_conv:
	fcvtns w2, d1
	fcvtnu x2, s1
	scvtf d1, w2
	ucvtf s1, x2
	fcvtps x2, d1
	fcvtms x2, d1
	fcvtzs x2, d1
	fcvtas x2, d1
	fcvtau x2, d1
	fjcvtzs w2, d1

fp_one_source:
	fmov h1, h2
	fmov s1, s2
	fmov d1, d2
	fabs d1, d2
	fneg d1, d2
	fsqrt d1, d2
	fcvt d1, s1
	fcvt d1, h1
	fcvt s1, d2
	frintn s1, s2
	frintp d1, d2
	frintm s1, s2
	frintz d1, d2
	frinta s1, s2
	frintx d1, d2
	frinti s1, s2
	frint32z d1, d2
	frint32x d1, d2
	frint64z d1, d2
	frint64x d1, d2

fp_compare:
	fcmp h1, h2
	fcmp s1, s2
	fcmp d1, d2
	fcmp s1, #0.0
	fcmp d1, #0.0
	fcmpe d1, d2
	fcmpe d1, #0.0

fp_immediate:
	fmov h1, #0.5
	fmov d1, #1.0
	fmov d1, #3.5
	fmov d1, #10.0
	fmov d1, #-1.5
	fmov d1, #21.0

fp_cond_comp:
	fccmp h1, h2, #0xe, eq
	fccmpe d1, d2, #0xa, ne

fp_two_source:
	fmul h0, h1, h2
	fmul s0, s1, s2
	fdiv d0, d1, d2
	fadd s0, s1, s2
	fsub d0, d1, d2
	fmax s0, s1, s2
	fmin d0, d1, d2
	fmaxnm s0, s1, s2
	fminnm d0, d1, d2
	fnmul s0, s1, s2

fp_cond_sel:
	fcsel h0, h1, h2, eq
	fcsel s0, s1, s2, ge
	fcsel d0, d1, d2, lt

fp_three_source:
	fmadd h0, h1, h2, h3
	fmsub s0, s1, s2, s3
	fnmadd d0, d1, d2, d3
	fnmsub s0, s1, s2, s3
