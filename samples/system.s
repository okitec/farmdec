exceptions:
	svc #1337
	hvc #3456
	smc #1234
	brk #42
	hlt #65535
	dcps1
	dcps2
	dcps3

hints:
	nop
	yield
	wfe
	wfi
	sev
	sevl
	esb
	tsb csync
	csdb
	bti

barriers:
	clrex
	dmb sy
	dmb st
	dmb ld
	dmb ish
	dmb ishst
	dmb ishld
	dmb nsh
	dmb nshst
	dmb nshld
	dmb osh
	dmb oshst
	dmb oshld
	isb sy
	sb
	dsb sy
	ssbb
	pssbb

pstate_inst:
	msr uao, #1
	msr pan, #1
	msr spsel, #0
	msr ssbs, #1
	msr dit, #1
	msr daifset, #5
	msr daifclr, #3
	cfinv
	xaflag
	axflag
	
sys_inst:
	at s1e1r, x23
	at s1e1w, x23
	at s1e0r, x23
	at s1e0w, x23
	at s1e1rp, x23
	at s1e1wp, x23
	at s1e2r, x23
	at s1e2w, x23
	at s12e1r, x23
	at s12e1w, x23
	at s12e0r, x23
	at s12e0w, x23
	at s1e3r, x23
	at s1e3w, x23
	cfp rctx, x12
	cpp rctx, x12
	dc ivac, x23
	dc isw, x23
	dc csw, x23
	dc cisw, x23
	dc zva, x23
	dc cvac, x23
	dc cvau, x23
	dc civac, x23
	dvp rctx, x12
	ic ialluis
	ic iallu
	ic ivau, x12
	tlbi vae3, x12 // I'm not enumerating two pages
	sysl x2, #1, c7, c8, #0

sys_reg_mov:
	msr nzcv, x2
	mrs x2, nzcv
