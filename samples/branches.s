cond_branch:
	b.eq cond_branch
	b.ne cond_branch
	b.cs cond_branch
	b.hs cond_branch
	b.cc cond_branch
	b.lo cond_branch
	b.mi cond_branch
	b.pl cond_branch
	b.vs cond_branch
	b.vc cond_branch
	b.hi cond_branch
	b.ls cond_branch
	b.ge cond_branch
	b.lt cond_branch
	b.gt cond_branch
	b.le cond_branch
	b.al cond_branch
	b.nv cond_branch

uncond_reg:
	br x3
	blr x3
	ret
	ret x3

uncond_imm:
	b  uncond_imm
	bl uncond_imm

cmp_and_branch:
	cbz x3, cmp_and_branch
	cbnz x3, cmp_and_branch

test_and_branch:
	tbz x3, #15, test_and_branch
	tbnz w1, #3, test_and_branch
