SIMD_mult:
	ld1 {v0.1d}, [sp]
	st1 {v0.2d}, [sp]
	ld1 {v0.16b, v1.16b, v2.16b}, [x0]
	st3 {v0.8b, v1.8b, v2.8b}, [sp]
	ld4 {v0.4s, v1.4s, v2.4s, v3.4s}, [sp], x3
	st2 {v0.2s, v1.2s}, [sp], #16

SIMD_single:
	st3 {v0.s, v1.s, v2.s}[3], [sp]
	ld1 {v0.b}[15], [sp], x3
	st4 {v0.d, v1.d, v2.d, v3.d}[0], [sp], #32
	ld2r {v0.8h, v1.8h}, [sp]
	ld4r {v0.4s, v1.4s, v2.4s, v3.4s}, [sp], x2

exclusive:
	ldxr x0, [sp]
	stxr w4, w1, [sp]
	stlxr  w4, w1, [sp]
	ldxrb w2, [sp]
	ldaxrb w2, [sp]
	ldxp x1, x2, [x23]
	ldaxp x1, x2, [x23]
	stxp w3, x4, x2, [x24]
	stlxp w3, x4, x2, [x24]
	casp w2, w3, w4, w5, [sp]
	casal x2, x5, [sp]

ordered:
	ldar w0, [sp]
	stlr w4, [sp]
	ldlar x2, [sp]
	stllr x3, [sp]
	ldapr x1, [sp]
	ldapur x2, [sp, #-3] // XXX "Tag instructions not supported"

no_alloc_pair:
	ldnp w2, w2, [sp, #-256]
	stnp x1, x2, [sp]
	ldnp s1, s2, [sp, #4]
	stnp d1, d2, [sp, #16]
	ldnp q1, q2, [sp, #-16]

pair:
	ldp w2, w2, [sp], #-256
	ldpsw x3, x5, [sp]
	stp x1, x2, [sp]
	ldp x1, x2, [sp, #16]!
	ldp s1, s2, [sp, #4]
	stp d1, d2, [sp, #16]
	ldp q1, q2, [sp, #-16]

unscaled:
	sturb w1, [x0, #-10]
	ldurb w2, [sp, #7]
	ldursb w3, [sp, #34]
	ldursb x0, [sp, #-19]
	stur w1, [x0, #-10]
	ldur w2, [sp, #7]
	ldur x3, [sp, #34]
	ldursw x0, [sp, #-19]
	ldur q0, [sp, #22]
	stur q1, [sp, #1]

register:
	ldr w4, register
	ldrsw x2, register
	ldr q5, register
	ldr w2, [sp], #-256
	str x1, [sp]
	ldr x1, [sp, #16]!
	ldr s1, [sp, #4]
	str d1, [sp, #16]
	ldr q1, [sp, #32]
	ldr x3, [sp, x4, lsl #3]
	str x2, [sp, w2, uxtw #0]

prefetch:
	prfm pldL1keep, [sp, #8]
	prfm pliL2strm, prefetch
	prfm pstL3keep, [sp, w0, sxtw #3]
	prfum pldL1strm, [sp, #-43]

// Just the LD* variants, since the ST* variants are just aliases.
atomic:
	ldaddb   w3, w5, [sp]
	ldaddab  w3, w5, [sp]
	ldaddlb  w3, w5, [sp]
	ldaddalb w3, w5, [sp]
	ldclr x2, x7, [x0]
	ldeor w5, w7, [x4]
	ldseth w3, w4, [sp]
	ldsmax x2, x6, [sp]
	ldsmin x2, x6, [sp]
	ldumax x2, x6, [sp]
	ldumin x2, x6, [sp]
	swp w2, w5, [x14]
	swp x3, x6, [x4]
