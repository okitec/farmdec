// maxArray example from Godbolt, compiled with -O0

maxArray:                               // @maxArray
	sub     sp, sp, #32             // =32
	str     x0, [sp, #24]
	str     x1, [sp, #16]
	str     wzr, [sp, #12]
.LBB0_1:                                // =>This Inner Loop Header: Depth=1
	ldr     w8, [sp, #12]
	cmp     w8, #16, lsl #12        // =65536
	b.ge    .LBB0_6
	ldr     x8, [sp, #16]
	ldrsw   x9, [sp, #12]
	ldr     d0, [x8, x9, lsl #3]
	ldr     x8, [sp, #24]
	ldrsw   x9, [sp, #12]
	ldr     d1, [x8, x9, lsl #3]
	fcmp    d0, d1
	b.le    .LBB0_4
	ldr     x8, [sp, #16]
	ldrsw   x9, [sp, #12]
	ldr     x8, [x8, x9, lsl #3]
	ldr     x9, [sp, #24]
	ldrsw   x10, [sp, #12]
	str     x8, [x9, x10, lsl #3]
.LBB0_4:                                //   in Loop: Header=BB0_1 Depth=1
	ldr     w8, [sp, #12]
	add     w8, w8, #1              // =1
	str     w8, [sp, #12]
	b       .LBB0_1
.LBB0_6:
	add     sp, sp, #32             // =32
	ret
