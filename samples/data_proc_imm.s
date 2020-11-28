// Examples of immediate data processing instructions.

pc_rel_addr:
	adr x0, extends
	adrp x0, imm_shifts

addsub:
	add sp, x4, #48
	adds x0, x1, #100, lsl #12
	cmn x3, #42
	sub w0, w1, #20
	subs x23, x24, #123
	cmp x3, #16

logical:
	and x0, x1, #0xff
	ands x0, x1, #0xff00
	tst x1, #0xff
	orr w0, w1, #0xff0000ff
	eor x0, x1, #0xf0f0f0f0f0f0f0f0

move_wide:
	mov x0, #1234
	movk w1, #1234
	movn x0, #0xaaaa
	movn x0, #0xaaaa, lsl #32
	movz x4, #0xff00
	movz x4, #0xff00, lsl #48

imm_shifts:
	lsl x0, x1, #32
	lsr x0, x1, #7
	asr x0, x1, #3
	ror x0, x1, #15

bitfield:
	ubfx x0, x1, #12, #3
	sbfx w0, w1, #7, #3
	ubfiz w0, w1, #4, #2
	sbfiz x0, x1, #15, #7
	bfc x0, #4, #3
	bfi x0, x1, #4, #3
	bfxil w0, w1, #23, #3

extends:
	uxtb w0, w1
	uxth w0, w1
	sxtb w0, w1
	sxth w0, w1
	sxtw x0, w1

extract:
	extr x0, x1, x2, #7

