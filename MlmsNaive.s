	.file	"MlmsNaive.cpp"
	.text
	.p2align 4
	.type	_Z16calculation_loopR6matrixRKS_ddd._omp_fn.0, @function
_Z16calculation_loopR6matrixRKS_ddd._omp_fn.0:
.LFB4426:
	.cfi_startproc
	endbr64
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$200, %rsp
	.cfi_def_cfa_offset 256
	movq	(%rdi), %rax
	movq	%rdi, 176(%rsp)
	movq	%rax, 152(%rsp)
	movq	(%rax), %rax
	movq	%rax, 120(%rsp)
	testq	%rax, %rax
	je	.L67
	movq	%rax, %r15
	call	omp_get_num_threads@PLT
	movl	%eax, %ebx
	call	omp_get_thread_num@PLT
	xorl	%edx, %edx
	movslq	%eax, %rcx
	movslq	%ebx, %rsi
	movq	%r15, %rax
	divq	%rsi
	cmpq	%rdx, %rcx
	jb	.L3
.L43:
	imulq	%rax, %rcx
	movq	%rax, %rdi
	addq	%rdx, %rcx
	addq	%rcx, %rdi
	movq	%rcx, 160(%rsp)
	movq	%rcx, %rax
	movq	%rdi, 184(%rsp)
	cmpq	%rdi, %rcx
	jnb	.L67
	movq	176(%rsp), %rcx
	vmovsd	32(%rcx), %xmm7
	vmovsd	16(%rcx), %xmm2
	vmovsd	%xmm7, 104(%rsp)
	vmovsd	24(%rcx), %xmm7
	movq	152(%rsp), %rcx
	vmovsd	%xmm2, 88(%rsp)
	movq	8(%rcx), %rbx
	vmovsd	%xmm7, 112(%rsp)
	testq	%rbx, %rbx
	je	.L67
	vmulsd	.LC1(%rip), %xmm2, %xmm1
	imulq	%rbx, %rax
	movq	%rax, 168(%rsp)
	vmovsd	%xmm1, 24(%rsp)
.L8:
	movq	176(%rsp), %rax
	movq	8(%rax), %r12
	movq	168(%rsp), %rax
	salq	$3, %rax
	movq	%rax, 144(%rsp)
	movq	160(%rsp), %rax
	testq	%rax, %rax
	js	.L5
	vxorpd	%xmm7, %xmm7, %xmm7
	vcvtsi2sdq	%rax, %xmm7, %xmm0
.L6:
	vmulsd	88(%rsp), %xmm0, %xmm2
	movq	$0, 136(%rsp)
	vmovsd	%xmm2, 128(%rsp)
	.p2align 4
	.p2align 3
.L7:
	movq	136(%rsp), %rax
	testq	%rax, %rax
	js	.L9
	vxorpd	%xmm1, %xmm1, %xmm1
	vcvtsi2sdq	%rax, %xmm1, %xmm0
.L10:
	vmulsd	88(%rsp), %xmm0, %xmm7
	xorl	%r13d, %r13d
	movq	%rbx, %rax
	xorl	%r15d, %r15d
	movq	%r13, %rbx
	movq	%rax, %r13
	vmovsd	%xmm7, 96(%rsp)
	.p2align 4
	.p2align 3
.L11:
	testq	%r15, %r15
	js	.L40
	vxorpd	%xmm1, %xmm1, %xmm1
	vcvtsi2sdq	%r15, %xmm1, %xmm0
.L41:
	vmovsd	128(%rsp), %xmm2
	vmovsd	24(%rsp), %xmm4
	xorl	%ebp, %ebp
	vfnmadd132sd	88(%rsp), %xmm2, %xmm0
	vaddsd	%xmm4, %xmm0, %xmm2
	vsubsd	%xmm4, %xmm0, %xmm3
	vmulsd	%xmm2, %xmm2, %xmm6
	vmovsd	%xmm2, 16(%rsp)
	vmovsd	%xmm3, (%rsp)
	vmovsd	%xmm6, 32(%rsp)
	jmp	.L38
	.p2align 4
	.p2align 3
.L71:
	vxorpd	%xmm7, %xmm7, %xmm7
	vcvtsi2sdq	%rbp, %xmm7, %xmm1
.L13:
	vmovsd	96(%rsp), %xmm7
	vfnmadd132sd	88(%rsp), %xmm7, %xmm1
	vxorpd	%xmm7, %xmm7, %xmm7
	vaddsd	24(%rsp), %xmm1, %xmm6
	vmulsd	%xmm6, %xmm6, %xmm4
	vaddsd	32(%rsp), %xmm4, %xmm2
	vmovsd	%xmm6, 8(%rsp)
	vucomisd	%xmm2, %xmm7
	ja	.L58
	vsqrtsd	%xmm2, %xmm2, %xmm3
.L16:
	vmovsd	(%rsp), %xmm6
	vaddsd	16(%rsp), %xmm3, %xmm3
	vmulsd	%xmm6, %xmm6, %xmm5
	vmovq	%xmm5, %r14
	vaddsd	%xmm4, %xmm5, %xmm5
	vmovsd	%xmm5, 40(%rsp)
	vucomisd	%xmm5, %xmm7
	ja	.L59
	vsqrtsd	%xmm5, %xmm5, %xmm0
.L19:
	vaddsd	(%rsp), %xmm0, %xmm0
	vmovsd	%xmm2, 64(%rsp)
	vmovsd	%xmm1, 56(%rsp)
	vdivsd	%xmm0, %xmm3, %xmm0
	call	log@PLT
	vmulsd	8(%rsp), %xmm0, %xmm3
	vxorpd	%xmm6, %xmm6, %xmm6
	vmovsd	64(%rsp), %xmm2
	vmovsd	56(%rsp), %xmm1
	vucomisd	%xmm2, %xmm6
	vmovsd	%xmm3, 48(%rsp)
	ja	.L60
	vsqrtsd	%xmm2, %xmm2, %xmm2
.L22:
	vsubsd	24(%rsp), %xmm1, %xmm1
	vaddsd	8(%rsp), %xmm2, %xmm2
	vmulsd	%xmm1, %xmm1, %xmm3
	vaddsd	32(%rsp), %xmm3, %xmm4
	vucomisd	%xmm4, %xmm6
	ja	.L61
	vsqrtsd	%xmm4, %xmm4, %xmm0
.L25:
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm4, 80(%rsp)
	vmovsd	%xmm3, 72(%rsp)
	vmovsd	%xmm1, 64(%rsp)
	vdivsd	%xmm0, %xmm2, %xmm0
	call	log@PLT
	vmovsd	72(%rsp), %xmm3
	vmovq	%r14, %xmm5
	vxorpd	%xmm7, %xmm7, %xmm7
	vmulsd	16(%rsp), %xmm0, %xmm4
	vmovsd	64(%rsp), %xmm1
	vaddsd	%xmm5, %xmm3, %xmm2
	vmovsd	%xmm4, 56(%rsp)
	vucomisd	%xmm2, %xmm7
	vmovsd	80(%rsp), %xmm4
	ja	.L62
	vmovsd	%xmm7, %xmm7, %xmm5
	vsqrtsd	%xmm2, %xmm2, %xmm3
.L28:
	vaddsd	(%rsp), %xmm3, %xmm3
	vucomisd	%xmm4, %xmm5
	ja	.L63
	vsqrtsd	%xmm4, %xmm4, %xmm4
.L31:
	vaddsd	16(%rsp), %xmm4, %xmm4
	vmovsd	%xmm2, 72(%rsp)
	vmovsd	%xmm1, 64(%rsp)
	vdivsd	%xmm4, %xmm3, %xmm0
	call	log@PLT
	vmovsd	64(%rsp), %xmm1
	vxorpd	%xmm3, %xmm3, %xmm3
	vmovsd	72(%rsp), %xmm2
	vmulsd	%xmm0, %xmm1, %xmm4
	vucomisd	%xmm2, %xmm3
	vmovq	%xmm4, %r14
	ja	.L64
	vmovsd	%xmm3, %xmm3, %xmm4
	vsqrtsd	%xmm2, %xmm2, %xmm2
.L34:
	vmovsd	40(%rsp), %xmm7
	vaddsd	%xmm2, %xmm1, %xmm1
	vucomisd	%xmm7, %xmm4
	ja	.L65
	vsqrtsd	%xmm7, %xmm7, %xmm0
.L37:
	vaddsd	8(%rsp), %xmm0, %xmm0
	vdivsd	%xmm0, %xmm1, %xmm0
	call	log@PLT
	vmovsd	48(%rsp), %xmm3
	vmovq	%r14, %xmm7
	vmovsd	%xmm0, %xmm0, %xmm1
	vaddsd	56(%rsp), %xmm3, %xmm0
	vmovsd	.LC2(%rip), %xmm2
	vmovsd	.LC3(%rip), %xmm3
	movq	8(%r12), %rax
	movq	16(%r12), %rdx
	imulq	%r15, %rax
	addq	%rbp, %rax
	incq	%rbp
	vaddsd	%xmm7, %xmm0, %xmm0
	vfmadd132sd	(%rsp), %xmm0, %xmm1
	vsubsd	112(%rsp), %xmm2, %xmm0
	vmovq	%rbx, %xmm7
	vmulsd	104(%rsp), %xmm3, %xmm2
	vdivsd	%xmm2, %xmm0, %xmm0
	vmulsd	%xmm1, %xmm0, %xmm0
	vfmadd231sd	(%rdx,%rax,8), %xmm0, %xmm7
	vmovq	%xmm7, %rbx
	cmpq	%r13, %rbp
	je	.L70
.L38:
	testq	%rbp, %rbp
	jns	.L71
	movq	%rbp, %rax
	movq	%rbp, %rdx
	vxorpd	%xmm4, %xmm4, %xmm4
	shrq	%rax
	andl	$1, %edx
	orq	%rdx, %rax
	vcvtsi2sdq	%rax, %xmm4, %xmm1
	vaddsd	%xmm1, %xmm1, %xmm1
	jmp	.L13
	.p2align 4
	.p2align 3
.L70:
	incq	%r15
	cmpq	%r15, 120(%rsp)
	jne	.L11
	movq	152(%rsp), %rax
	movq	144(%rsp), %rsi
	movq	%r13, %rbx
	incq	136(%rsp)
	movq	16(%rax), %rax
	vmovsd	%xmm7, (%rax,%rsi)
	movq	136(%rsp), %rax
	addq	$8, %rsi
	movq	%rsi, 144(%rsp)
	cmpq	%r13, %rax
	jne	.L7
	incq	160(%rsp)
	addq	%rbx, 168(%rsp)
	movq	160(%rsp), %rax
	cmpq	%rax, 184(%rsp)
	jne	.L8
.L67:
	addq	$200, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4
	.p2align 3
.L40:
	.cfi_restore_state
	movq	%r15, %rax
	movq	%r15, %rdx
	vxorpd	%xmm2, %xmm2, %xmm2
	shrq	%rax
	andl	$1, %edx
	orq	%rdx, %rax
	vcvtsi2sdq	%rax, %xmm2, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm0
	jmp	.L41
.L9:
	movq	%rax, %rdi
	shrq	%rax
	vxorpd	%xmm2, %xmm2, %xmm2
	andl	$1, %edi
	orq	%rdi, %rax
	vcvtsi2sdq	%rax, %xmm2, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm0
	jmp	.L10
.L5:
	movq	%rax, %rcx
	shrq	%rax
	vxorpd	%xmm1, %xmm1, %xmm1
	andl	$1, %ecx
	orq	%rcx, %rax
	vcvtsi2sdq	%rax, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm0
	jmp	.L6
.L65:
	vmovsd	%xmm7, %xmm7, %xmm0
	vmovsd	%xmm1, 64(%rsp)
	call	sqrt@PLT
	vmovsd	64(%rsp), %xmm1
	jmp	.L37
.L64:
	vmovsd	%xmm2, %xmm2, %xmm0
	vmovsd	%xmm1, 64(%rsp)
	call	sqrt@PLT
	vxorpd	%xmm4, %xmm4, %xmm4
	vmovsd	64(%rsp), %xmm1
	vmovsd	%xmm0, %xmm0, %xmm2
	jmp	.L34
.L63:
	vmovsd	%xmm4, %xmm4, %xmm0
	vmovsd	%xmm3, 80(%rsp)
	vmovsd	%xmm2, 72(%rsp)
	vmovsd	%xmm1, 64(%rsp)
	call	sqrt@PLT
	vmovsd	80(%rsp), %xmm3
	vmovsd	72(%rsp), %xmm2
	vmovsd	64(%rsp), %xmm1
	vmovsd	%xmm0, %xmm0, %xmm4
	jmp	.L31
.L62:
	vmovsd	%xmm2, %xmm2, %xmm0
	vmovsd	%xmm4, 80(%rsp)
	vmovsd	%xmm1, 72(%rsp)
	vmovsd	%xmm2, 64(%rsp)
	call	sqrt@PLT
	vmovsd	80(%rsp), %xmm4
	vxorpd	%xmm5, %xmm5, %xmm5
	vmovsd	72(%rsp), %xmm1
	vmovsd	64(%rsp), %xmm2
	vmovsd	%xmm0, %xmm0, %xmm3
	jmp	.L28
.L61:
	vmovsd	%xmm4, %xmm4, %xmm0
	vmovsd	%xmm1, 80(%rsp)
	vmovsd	%xmm2, 72(%rsp)
	vmovsd	%xmm3, 64(%rsp)
	vmovsd	%xmm4, 56(%rsp)
	call	sqrt@PLT
	vmovsd	80(%rsp), %xmm1
	vmovsd	72(%rsp), %xmm2
	vmovsd	64(%rsp), %xmm3
	vmovsd	56(%rsp), %xmm4
	jmp	.L25
.L60:
	vmovsd	%xmm2, %xmm2, %xmm0
	vmovsd	%xmm1, 56(%rsp)
	call	sqrt@PLT
	vxorpd	%xmm6, %xmm6, %xmm6
	vmovsd	56(%rsp), %xmm1
	vmovsd	%xmm0, %xmm0, %xmm2
	jmp	.L22
.L59:
	vmovsd	%xmm5, %xmm5, %xmm0
	vmovsd	%xmm3, 64(%rsp)
	vmovsd	%xmm2, 56(%rsp)
	vmovsd	%xmm1, 48(%rsp)
	call	sqrt@PLT
	vmovsd	64(%rsp), %xmm3
	vmovsd	56(%rsp), %xmm2
	vmovsd	48(%rsp), %xmm1
	jmp	.L19
.L58:
	vmovsd	%xmm2, %xmm2, %xmm0
	vmovsd	%xmm4, 56(%rsp)
	vmovsd	%xmm1, 48(%rsp)
	vmovsd	%xmm2, 40(%rsp)
	call	sqrt@PLT
	vmovsd	56(%rsp), %xmm4
	vxorpd	%xmm7, %xmm7, %xmm7
	vmovsd	48(%rsp), %xmm1
	vmovsd	40(%rsp), %xmm2
	vmovsd	%xmm0, %xmm0, %xmm3
	jmp	.L16
.L3:
	incq	%rax
	xorl	%edx, %edx
	jmp	.L43
	.cfi_endproc
.LFE4426:
	.size	_Z16calculation_loopR6matrixRKS_ddd._omp_fn.0, .-_Z16calculation_loopR6matrixRKS_ddd._omp_fn.0
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC4:
	.string	"%.1f\t"
.LC5:
	.string	"\n"
	.text
	.p2align 4
	.globl	_Z10printarrayRK6matrix
	.type	_Z10printarrayRK6matrix, @function
_Z10printarrayRK6matrix:
.LFB3635:
	.cfi_startproc
	endbr64
	cmpq	$0, (%rdi)
	je	.L87
	pushq	%r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	movq	%rdi, %r14
	pushq	%r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	pushq	%r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	leaq	.LC5(%rip), %r12
	pushq	%rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	leaq	.LC4(%rip), %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L73:
	movq	8(%r14), %rax
	xorl	%r13d, %r13d
	testq	%rax, %rax
	je	.L76
	.p2align 4
	.p2align 3
.L75:
	movq	16(%r14), %rdx
	imulq	%rbx, %rax
	movq	%rbp, %rsi
	movl	$1, %edi
	addq	%r13, %rax
	incq	%r13
	vmovsd	(%rdx,%rax,8), %xmm0
	movl	$1, %eax
	call	__printf_chk@PLT
	movq	8(%r14), %rax
	cmpq	%r13, %rax
	ja	.L75
.L76:
	movq	%r12, %rsi
	movl	$1, %edi
	xorl	%eax, %eax
	incq	%rbx
	call	__printf_chk@PLT
	cmpq	%rbx, (%r14)
	ja	.L73
	popq	%rbx
	.cfi_def_cfa_offset 40
	popq	%rbp
	.cfi_def_cfa_offset 32
	popq	%r12
	.cfi_def_cfa_offset 24
	popq	%r13
	.cfi_def_cfa_offset 16
	popq	%r14
	.cfi_def_cfa_offset 8
	ret
.L87:
	.cfi_restore 3
	.cfi_restore 6
	.cfi_restore 12
	.cfi_restore 13
	.cfi_restore 14
	ret
	.cfi_endproc
.LFE3635:
	.size	_Z10printarrayRK6matrix, .-_Z10printarrayRK6matrix
	.p2align 4
	.globl	_Z23initializePressureArrayR6matrixddd
	.type	_Z23initializePressureArrayR6matrixddd, @function
_Z23initializePressureArrayR6matrixddd:
.LFB3638:
	.cfi_startproc
	endbr64
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	vmovq	%xmm0, %r13
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rdi, %rbp
	subq	$24, %rsp
	.cfi_def_cfa_offset 80
	movq	(%rdi), %r12
	testq	%r12, %r12
	je	.L89
	movq	8(%rdi), %r15
	testq	%r15, %r15
	je	.L89
	salq	$3, %r15
	xorl	%r14d, %r14d
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L90:
	movq	16(%rbp), %rdi
	xorl	%esi, %esi
	movq	%r15, %rdx
	incq	%rbx
	vmovsd	%xmm2, 8(%rsp)
	vmovsd	%xmm1, (%rsp)
	addq	%r14, %rdi
	addq	%r15, %r14
	call	memset@PLT
	cmpq	%r12, %rbx
	vmovsd	(%rsp), %xmm1
	vmovsd	8(%rsp), %xmm2
	jne	.L90
.L89:
	vmovsd	.LC6(%rip), %xmm0
	vmovq	%r13, %xmm4
	vcomisd	%xmm0, %xmm4
	jnb	.L91
	vxorps	%xmm3, %xmm3, %xmm3
	vcvttsd2siq	%xmm4, %r8
	testq	%r8, %r8
	js	.L94
.L117:
	vcvtsi2sdq	%r8, %xmm3, %xmm0
.L95:
	vcomisd	%xmm0, %xmm1
	jbe	.L113
	movq	8(%rbp), %r10
	movq	16(%rbp), %r9
	movq	%r8, %rdi
	.p2align 4
	.p2align 3
.L97:
	movq	%r10, %rax
	imulq	%rdi, %rax
	leaq	(%r9,%rax,8), %rcx
	movq	%r8, %rax
	jmp	.L98
	.p2align 4
	.p2align 3
.L116:
	vcvtsi2sdq	%rax, %xmm3, %xmm0
	vcomisd	%xmm0, %xmm1
	jbe	.L115
.L98:
	vmovsd	%xmm2, (%rcx,%rax,8)
	incq	%rax
	jns	.L116
	movq	%rax, %rdx
	movq	%rax, %rsi
	shrq	%rdx
	andl	$1, %esi
	orq	%rsi, %rdx
	vcvtsi2sdq	%rdx, %xmm3, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm0
	vcomisd	%xmm0, %xmm1
	ja	.L98
.L115:
	incq	%rdi
	js	.L101
	vcvtsi2sdq	%rdi, %xmm3, %xmm0
.L102:
	vcomisd	%xmm0, %xmm1
	ja	.L97
.L113:
	addq	$24, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4
	.p2align 3
.L101:
	.cfi_restore_state
	movq	%rdi, %rax
	movq	%rdi, %rdx
	shrq	%rax
	andl	$1, %edx
	orq	%rdx, %rax
	vcvtsi2sdq	%rax, %xmm3, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm0
	jmp	.L102
.L91:
	vsubsd	%xmm0, %xmm4, %xmm0
	vxorps	%xmm3, %xmm3, %xmm3
	vcvttsd2siq	%xmm0, %r8
	btcq	$63, %r8
	testq	%r8, %r8
	jns	.L117
.L94:
	movq	%r8, %rax
	movq	%r8, %rdx
	shrq	%rax
	andl	$1, %edx
	orq	%rdx, %rax
	vcvtsi2sdq	%rax, %xmm3, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm0
	jmp	.L95
	.cfi_endproc
.LFE3638:
	.size	_Z23initializePressureArrayR6matrixddd, .-_Z23initializePressureArrayR6matrixddd
	.p2align 4
	.globl	_Z27initializeDisplacementArrayR6matrix
	.type	_Z27initializeDisplacementArrayR6matrix, @function
_Z27initializeDisplacementArrayR6matrix:
.LFB3639:
	.cfi_startproc
	endbr64
	pushq	%r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	pushq	%r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	pushq	%r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	pushq	%rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	pushq	%rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	movq	(%rdi), %r14
	testq	%r14, %r14
	je	.L128
	movq	8(%rdi), %r12
	movq	%rdi, %r13
	testq	%r12, %r12
	je	.L128
	salq	$3, %r12
	xorl	%ebp, %ebp
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L120:
	movq	16(%r13), %rdi
	movq	%r12, %rdx
	xorl	%esi, %esi
	incq	%rbx
	addq	%rbp, %rdi
	addq	%r12, %rbp
	call	memset@PLT
	cmpq	%r14, %rbx
	jne	.L120
.L128:
	popq	%rbx
	.cfi_def_cfa_offset 40
	popq	%rbp
	.cfi_def_cfa_offset 32
	popq	%r12
	.cfi_def_cfa_offset 24
	popq	%r13
	.cfi_def_cfa_offset 16
	popq	%r14
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE3639:
	.size	_Z27initializeDisplacementArrayR6matrix, .-_Z27initializeDisplacementArrayR6matrix
	.p2align 4
	.globl	_Z16calculation_loopR6matrixRKS_ddd
	.type	_Z16calculation_loopR6matrixRKS_ddd, @function
_Z16calculation_loopR6matrixRKS_ddd:
.LFB3642:
	.cfi_startproc
	endbr64
	subq	$56, %rsp
	.cfi_def_cfa_offset 64
	vmovq	%rdi, %xmm3
	vunpcklpd	%xmm1, %xmm0, %xmm1
	xorl	%ecx, %ecx
	vpinsrq	$1, %rsi, %xmm3, %xmm0
	xorl	%edx, %edx
	movq	%rsp, %rsi
	leaq	_Z16calculation_loopR6matrixRKS_ddd._omp_fn.0(%rip), %rdi
	movq	%fs:40, %rax
	movq	%rax, 40(%rsp)
	xorl	%eax, %eax
	vmovsd	%xmm2, 32(%rsp)
	vmovapd	%xmm1, 16(%rsp)
	vmovdqa	%xmm0, (%rsp)
	call	GOMP_parallel@PLT
	movq	40(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L134
	addq	$56, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
.L134:
	.cfi_restore_state
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE3642:
	.size	_Z16calculation_loopR6matrixRKS_ddd, .-_Z16calculation_loopR6matrixRKS_ddd
	.p2align 4
	.globl	_Z17calc_displacementRK6matrixS1_ddddddd
	.type	_Z17calc_displacementRK6matrixS1_ddddddd, @function
_Z17calc_displacementRK6matrixS1_ddddddd:
.LFB3643:
	.cfi_startproc
	endbr64
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$168, %rsp
	.cfi_def_cfa_offset 224
	movq	(%rsi), %rax
	vmovsd	%xmm0, 136(%rsp)
	vmovsd	%xmm1, 144(%rsp)
	vmovsd	%xmm2, 40(%rsp)
	vmovsd	%xmm3, 152(%rsp)
	vmovsd	%xmm4, 112(%rsp)
	vmovsd	%xmm5, 120(%rsp)
	vmovsd	%xmm6, 96(%rsp)
	movq	%rax, 128(%rsp)
	testq	%rax, %rax
	je	.L169
	movq	8(%rsi), %r15
	movq	%rdi, %r12
	xorl	%r14d, %r14d
	xorl	%r13d, %r13d
	.p2align 4
	.p2align 3
.L137:
	testq	%r15, %r15
	je	.L168
	vmovsd	144(%rsp), %xmm4
	vmulsd	96(%rsp), %xmm4, %xmm4
	vmovsd	%xmm4, 104(%rsp)
	testq	%r14, %r14
	js	.L138
	vxorpd	%xmm2, %xmm2, %xmm2
	vcvtsi2sdq	%r14, %xmm2, %xmm0
.L139:
	vmovsd	96(%rsp), %xmm2
	vmovsd	152(%rsp), %xmm3
	xorl	%ebp, %ebp
	vmulsd	%xmm2, %xmm0, %xmm0
	vfmsub231sd	136(%rsp), %xmm2, %xmm0
	vaddsd	%xmm3, %xmm0, %xmm1
	vsubsd	%xmm3, %xmm0, %xmm3
	vmulsd	%xmm1, %xmm1, %xmm4
	vmovsd	%xmm1, 24(%rsp)
	vmovsd	%xmm3, 8(%rsp)
	vmovsd	%xmm4, 32(%rsp)
	jmp	.L166
	.p2align 4
	.p2align 3
.L188:
	vxorpd	%xmm4, %xmm4, %xmm4
	vcvtsi2sdq	%rbp, %xmm4, %xmm1
.L141:
	vmovsd	104(%rsp), %xmm5
	vfnmadd132sd	96(%rsp), %xmm5, %xmm1
	vaddsd	40(%rsp), %xmm1, %xmm5
	vmulsd	%xmm5, %xmm5, %xmm4
	vaddsd	32(%rsp), %xmm4, %xmm2
	vmovsd	%xmm5, 16(%rsp)
	vxorpd	%xmm5, %xmm5, %xmm5
	vucomisd	%xmm2, %xmm5
	ja	.L179
	vsqrtsd	%xmm2, %xmm2, %xmm3
.L144:
	vmovsd	8(%rsp), %xmm7
	vaddsd	24(%rsp), %xmm3, %xmm3
	vmulsd	%xmm7, %xmm7, %xmm6
	vaddsd	%xmm6, %xmm4, %xmm7
	vmovsd	%xmm6, 64(%rsp)
	vmovsd	%xmm7, 48(%rsp)
	vucomisd	%xmm7, %xmm5
	ja	.L180
	vsqrtsd	%xmm7, %xmm7, %xmm0
.L147:
	vaddsd	8(%rsp), %xmm0, %xmm0
	vmovsd	%xmm2, 72(%rsp)
	vmovsd	%xmm1, 56(%rsp)
	vdivsd	%xmm0, %xmm3, %xmm0
	call	log@PLT
	vmulsd	16(%rsp), %xmm0, %xmm3
	vxorpd	%xmm7, %xmm7, %xmm7
	vmovsd	72(%rsp), %xmm2
	vmovsd	56(%rsp), %xmm1
	vucomisd	%xmm2, %xmm7
	vmovq	%xmm3, %rbx
	ja	.L181
	vsqrtsd	%xmm2, %xmm2, %xmm2
.L150:
	vsubsd	40(%rsp), %xmm1, %xmm1
	vaddsd	16(%rsp), %xmm2, %xmm2
	vmulsd	%xmm1, %xmm1, %xmm3
	vaddsd	32(%rsp), %xmm3, %xmm4
	vucomisd	%xmm4, %xmm7
	ja	.L182
	vsqrtsd	%xmm4, %xmm4, %xmm0
.L153:
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm4, 88(%rsp)
	vmovsd	%xmm3, 80(%rsp)
	vmovsd	%xmm1, 72(%rsp)
	vdivsd	%xmm0, %xmm2, %xmm0
	call	log@PLT
	vmulsd	24(%rsp), %xmm0, %xmm3
	vxorpd	%xmm4, %xmm4, %xmm4
	vmovsd	72(%rsp), %xmm1
	vmovsd	%xmm3, 56(%rsp)
	vmovsd	80(%rsp), %xmm3
	vaddsd	64(%rsp), %xmm3, %xmm2
	vucomisd	%xmm2, %xmm4
	vmovsd	88(%rsp), %xmm4
	ja	.L183
	vsqrtsd	%xmm2, %xmm2, %xmm3
.L156:
	vxorpd	%xmm6, %xmm6, %xmm6
	vaddsd	8(%rsp), %xmm3, %xmm3
	vucomisd	%xmm4, %xmm6
	ja	.L184
	vsqrtsd	%xmm4, %xmm4, %xmm4
.L159:
	vaddsd	24(%rsp), %xmm4, %xmm4
	vmovsd	%xmm2, 72(%rsp)
	vmovsd	%xmm1, 64(%rsp)
	vdivsd	%xmm4, %xmm3, %xmm0
	call	log@PLT
	vmovsd	64(%rsp), %xmm1
	vxorpd	%xmm4, %xmm4, %xmm4
	vmovsd	72(%rsp), %xmm2
	vmulsd	%xmm0, %xmm1, %xmm3
	vucomisd	%xmm2, %xmm4
	vmovsd	%xmm3, 64(%rsp)
	ja	.L185
	vmovsd	%xmm4, %xmm4, %xmm3
	vsqrtsd	%xmm2, %xmm2, %xmm2
.L162:
	vmovsd	48(%rsp), %xmm6
	vaddsd	%xmm2, %xmm1, %xmm1
	vucomisd	%xmm6, %xmm3
	ja	.L186
	vsqrtsd	%xmm6, %xmm6, %xmm0
.L165:
	vaddsd	16(%rsp), %xmm0, %xmm0
	vdivsd	%xmm0, %xmm1, %xmm0
	call	log@PLT
	vmovq	%rbx, %xmm6
	vmovsd	.LC2(%rip), %xmm3
	vmovq	%r13, %xmm5
	vmovsd	%xmm0, %xmm0, %xmm1
	vaddsd	56(%rsp), %xmm6, %xmm0
	vmovsd	.LC3(%rip), %xmm6
	vmulsd	120(%rsp), %xmm6, %xmm2
	movq	8(%r12), %rax
	movq	16(%r12), %rdx
	vaddsd	64(%rsp), %xmm0, %xmm0
	imulq	%r14, %rax
	vfmadd132sd	8(%rsp), %xmm0, %xmm1
	vsubsd	112(%rsp), %xmm3, %xmm0
	addq	%rbp, %rax
	incq	%rbp
	vdivsd	%xmm2, %xmm0, %xmm0
	vmulsd	%xmm1, %xmm0, %xmm0
	vfmadd231sd	(%rdx,%rax,8), %xmm0, %xmm5
	vmovq	%xmm5, %r13
	cmpq	%rbp, %r15
	je	.L168
.L166:
	testq	%rbp, %rbp
	jns	.L188
	movq	%rbp, %rax
	movq	%rbp, %rdx
	vxorpd	%xmm2, %xmm2, %xmm2
	shrq	%rax
	andl	$1, %edx
	orq	%rdx, %rax
	vcvtsi2sdq	%rax, %xmm2, %xmm1
	vaddsd	%xmm1, %xmm1, %xmm1
	jmp	.L141
	.p2align 4
	.p2align 3
.L168:
	incq	%r14
	cmpq	128(%rsp), %r14
	jne	.L137
.L135:
	addq	$168, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	vmovq	%r13, %xmm0
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4
	.p2align 3
.L138:
	.cfi_restore_state
	movq	%r14, %rax
	movq	%r14, %rdx
	vxorpd	%xmm1, %xmm1, %xmm1
	shrq	%rax
	andl	$1, %edx
	orq	%rdx, %rax
	vcvtsi2sdq	%rax, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm0
	jmp	.L139
.L169:
	movq	.LC0(%rip), %r13
	jmp	.L135
.L179:
	vmovsd	%xmm2, %xmm2, %xmm0
	vmovsd	%xmm4, 64(%rsp)
	vmovsd	%xmm1, 56(%rsp)
	vmovsd	%xmm2, 48(%rsp)
	call	sqrt@PLT
	vmovsd	64(%rsp), %xmm4
	vxorpd	%xmm5, %xmm5, %xmm5
	vmovsd	56(%rsp), %xmm1
	vmovsd	48(%rsp), %xmm2
	vmovsd	%xmm0, %xmm0, %xmm3
	jmp	.L144
.L186:
	vmovsd	%xmm6, %xmm6, %xmm0
	vmovsd	%xmm1, 72(%rsp)
	call	sqrt@PLT
	vmovsd	72(%rsp), %xmm1
	jmp	.L165
.L185:
	vmovsd	%xmm2, %xmm2, %xmm0
	vmovsd	%xmm1, 72(%rsp)
	call	sqrt@PLT
	vxorpd	%xmm3, %xmm3, %xmm3
	vmovsd	72(%rsp), %xmm1
	vmovsd	%xmm0, %xmm0, %xmm2
	jmp	.L162
.L184:
	vmovsd	%xmm4, %xmm4, %xmm0
	vmovsd	%xmm3, 80(%rsp)
	vmovsd	%xmm2, 72(%rsp)
	vmovsd	%xmm1, 64(%rsp)
	call	sqrt@PLT
	vmovsd	80(%rsp), %xmm3
	vmovsd	72(%rsp), %xmm2
	vmovsd	64(%rsp), %xmm1
	vmovsd	%xmm0, %xmm0, %xmm4
	jmp	.L159
.L183:
	vmovsd	%xmm2, %xmm2, %xmm0
	vmovsd	%xmm4, 80(%rsp)
	vmovsd	%xmm1, 72(%rsp)
	vmovsd	%xmm2, 64(%rsp)
	call	sqrt@PLT
	vmovsd	80(%rsp), %xmm4
	vmovsd	72(%rsp), %xmm1
	vmovsd	64(%rsp), %xmm2
	vmovsd	%xmm0, %xmm0, %xmm3
	jmp	.L156
.L182:
	vmovsd	%xmm4, %xmm4, %xmm0
	vmovsd	%xmm1, 88(%rsp)
	vmovsd	%xmm2, 80(%rsp)
	vmovsd	%xmm3, 72(%rsp)
	vmovsd	%xmm4, 56(%rsp)
	call	sqrt@PLT
	vmovsd	88(%rsp), %xmm1
	vmovsd	80(%rsp), %xmm2
	vmovsd	72(%rsp), %xmm3
	vmovsd	56(%rsp), %xmm4
	jmp	.L153
.L181:
	vmovsd	%xmm2, %xmm2, %xmm0
	vmovsd	%xmm1, 56(%rsp)
	call	sqrt@PLT
	vxorpd	%xmm7, %xmm7, %xmm7
	vmovsd	56(%rsp), %xmm1
	vmovsd	%xmm0, %xmm0, %xmm2
	jmp	.L150
.L180:
	vmovsd	%xmm7, %xmm7, %xmm0
	vmovsd	%xmm3, 80(%rsp)
	vmovsd	%xmm2, 72(%rsp)
	vmovsd	%xmm1, 56(%rsp)
	call	sqrt@PLT
	vmovsd	80(%rsp), %xmm3
	vmovsd	72(%rsp), %xmm2
	vmovsd	56(%rsp), %xmm1
	jmp	.L147
	.cfi_endproc
.LFE3643:
	.size	_Z17calc_displacementRK6matrixS1_ddddddd, .-_Z17calc_displacementRK6matrixS1_ddddddd
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.type	_GLOBAL__sub_I__Z10printarrayRK6matrix, @function
_GLOBAL__sub_I__Z10printarrayRK6matrix:
.LFB4425:
	.cfi_startproc
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	leaq	_ZStL8__ioinit(%rip), %rbp
	movq	%rbp, %rdi
	call	_ZNSt8ios_base4InitC1Ev@PLT
	movq	%rbp, %rsi
	movq	_ZNSt8ios_base4InitD1Ev@GOTPCREL(%rip), %rdi
	leaq	__dso_handle(%rip), %rdx
	popq	%rbp
	.cfi_def_cfa_offset 8
	jmp	__cxa_atexit@PLT
	.cfi_endproc
.LFE4425:
	.size	_GLOBAL__sub_I__Z10printarrayRK6matrix, .-_GLOBAL__sub_I__Z10printarrayRK6matrix
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I__Z10printarrayRK6matrix
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC0:
	.long	0
	.long	0
	.align 8
.LC1:
	.long	0
	.long	1071644672
	.align 8
.LC2:
	.long	0
	.long	1072693248
	.align 8
.LC3:
	.long	1413754136
	.long	1074340347
	.align 8
.LC6:
	.long	0
	.long	1138753536
	.hidden	__dso_handle
	.ident	"GCC: (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
