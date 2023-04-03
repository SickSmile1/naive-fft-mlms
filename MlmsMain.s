	.file	"MlmsMain.cpp"
	.text
	.section	.rodata._ZNSt6vectorIdSaIdEE17_M_default_appendEm.str1.1,"aMS",@progbits,1
.LC1:
	.string	"vector::_M_default_append"
	.section	.text._ZNSt6vectorIdSaIdEE17_M_default_appendEm,"axG",@progbits,_ZNSt6vectorIdSaIdEE17_M_default_appendEm,comdat
	.align 2
	.p2align 4
	.weak	_ZNSt6vectorIdSaIdEE17_M_default_appendEm
	.type	_ZNSt6vectorIdSaIdEE17_M_default_appendEm, @function
_ZNSt6vectorIdSaIdEE17_M_default_appendEm:
.LFB2038:
	.cfi_startproc
	endbr64
	testq	%rsi, %rsi
	je	.L31
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	movabsq	$1152921504606846975, %rax
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rdi, %r12
	subq	$24, %rsp
	.cfi_def_cfa_offset 80
	movq	8(%rdi), %rdx
	movq	(%rdi), %r14
	movq	%rsi, %rbx
	movq	%rdx, %rbp
	subq	%r14, %rbp
	movq	%rbp, %r13
	sarq	$3, %r13
	subq	%r13, %rax
	movq	%rax, %rcx
	movq	16(%rdi), %rax
	subq	%rdx, %rax
	sarq	$3, %rax
	cmpq	%rax, %rsi
	jbe	.L35
	cmpq	%rsi, %rcx
	jb	.L36
	cmpq	%r13, %rsi
	movq	%r13, %rax
	cmovnb	%rsi, %rax
	addq	%r13, %rax
	jc	.L7
	testq	%rax, %rax
	jne	.L37
	movq	%rbp, %r8
	xorl	%r15d, %r15d
	xorl	%ecx, %ecx
.L9:
	movq	%rbx, %rdx
	addq	%rcx, %rbp
	decq	%rdx
	movq	$0x000000000, 0(%rbp)
	je	.L13
	leaq	8(%rbp), %rdi
	salq	$3, %rdx
	xorl	%esi, %esi
	movq	%r8, 8(%rsp)
	movq	%rcx, (%rsp)
	call	memset@PLT
	movq	(%rsp), %rcx
	movq	8(%rsp), %r8
.L13:
	testq	%r8, %r8
	jg	.L38
	testq	%r14, %r14
	jne	.L39
.L15:
	addq	%r13, %rbx
	vmovq	%rcx, %xmm1
	movq	%r15, 16(%r12)
	leaq	(%rcx,%rbx,8), %rax
	vpinsrq	$1, %rax, %xmm1, %xmm0
	vmovdqu	%xmm0, (%r12)
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
.L35:
	.cfi_restore_state
	decq	%rbx
	movq	$0x000000000, (%rdx)
	leaq	8(%rdx), %rcx
	je	.L4
	leaq	(%rcx,%rbx,8), %rax
	movq	%rcx, %rdi
	xorl	%esi, %esi
	subq	%rdx, %rax
	leaq	-8(%rax), %rbx
	movq	%rbx, %rdx
	call	memset@PLT
	movq	%rax, %rcx
	addq	%rbx, %rcx
.L4:
	movq	%rcx, 8(%r12)
	addq	$24, %rsp
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
.L31:
	.cfi_restore 3
	.cfi_restore 6
	.cfi_restore 12
	.cfi_restore 13
	.cfi_restore 14
	.cfi_restore 15
	ret
	.p2align 4
	.p2align 3
.L38:
	.cfi_def_cfa_offset 80
	.cfi_offset 3, -56
	.cfi_offset 6, -48
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	.cfi_offset 14, -24
	.cfi_offset 15, -16
	movq	%r14, %rsi
	movq	%rcx, %rdi
	movq	%r8, %rdx
	call	memmove@PLT
	movq	16(%r12), %rsi
	movq	%rax, %rcx
	subq	%r14, %rsi
.L14:
	movq	%r14, %rdi
	movq	%rcx, (%rsp)
	call	_ZdlPvm@PLT
	movq	(%rsp), %rcx
	jmp	.L15
	.p2align 4
	.p2align 3
.L39:
	movq	16(%r12), %rsi
	subq	%r14, %rsi
	jmp	.L14
.L37:
	movabsq	$1152921504606846975, %rdx
	cmpq	%rdx, %rax
	cmova	%rdx, %rax
	leaq	0(,%rax,8), %r15
.L8:
	movq	%r15, %rdi
	call	_Znwm@PLT
	movq	(%r12), %r14
	movq	8(%r12), %r8
	movq	%rax, %rcx
	addq	%rax, %r15
	subq	%r14, %r8
	jmp	.L9
.L7:
	movabsq	$9223372036854775800, %r15
	jmp	.L8
.L36:
	leaq	.LC1(%rip), %rdi
	call	_ZSt20__throw_length_errorPKc@PLT
	.cfi_endproc
.LFE2038:
	.size	_ZNSt6vectorIdSaIdEE17_M_default_appendEm, .-_ZNSt6vectorIdSaIdEE17_M_default_appendEm
	.section	.text._ZN6matrixC2ESt5arrayImLm2EE,"axG",@progbits,_ZN6matrixC5ESt5arrayImLm2EE,comdat
	.align 2
	.p2align 4
	.weak	_ZN6matrixC2ESt5arrayImLm2EE
	.type	_ZN6matrixC2ESt5arrayImLm2EE, @function
_ZN6matrixC2ESt5arrayImLm2EE:
.LFB1353:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA1353
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	vpxor	%xmm0, %xmm0, %xmm0
	movq	%rdi, %rbx
	subq	$8, %rsp
	.cfi_def_cfa_offset 32
	movq	%rsi, (%rdi)
	movq	%rdx, 8(%rdi)
	vmovdqu	%xmm0, 16(%rdi)
	imulq	%rdx, %rsi
	movq	$0, 32(%rdi)
	testq	%rsi, %rsi
	jne	.L52
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.p2align 4
	.p2align 3
.L52:
	.cfi_restore_state
	addq	$16, %rdi
.LEHB0:
	call	_ZNSt6vectorIdSaIdEE17_M_default_appendEm
.LEHE0:
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
.L45:
	.cfi_restore_state
	endbr64
	movq	%rax, %rbp
.L43:
	movq	16(%rbx), %rdi
	movq	32(%rbx), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	je	.L49
	vzeroupper
	call	_ZdlPvm@PLT
.L44:
	movq	%rbp, %rdi
.LEHB1:
	call	_Unwind_Resume@PLT
.LEHE1:
.L49:
	vzeroupper
	jmp	.L44
	.cfi_endproc
.LFE1353:
	.globl	__gxx_personality_v0
	.section	.gcc_except_table._ZN6matrixC2ESt5arrayImLm2EE,"aG",@progbits,_ZN6matrixC5ESt5arrayImLm2EE,comdat
.LLSDA1353:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE1353-.LLSDACSB1353
.LLSDACSB1353:
	.uleb128 .LEHB0-.LFB1353
	.uleb128 .LEHE0-.LEHB0
	.uleb128 .L45-.LFB1353
	.uleb128 0
	.uleb128 .LEHB1-.LFB1353
	.uleb128 .LEHE1-.LEHB1
	.uleb128 0
	.uleb128 0
.LLSDACSE1353:
	.section	.text._ZN6matrixC2ESt5arrayImLm2EE,"axG",@progbits,_ZN6matrixC5ESt5arrayImLm2EE,comdat
	.size	_ZN6matrixC2ESt5arrayImLm2EE, .-_ZN6matrixC2ESt5arrayImLm2EE
	.weak	_ZN6matrixC1ESt5arrayImLm2EE
	.set	_ZN6matrixC1ESt5arrayImLm2EE,_ZN6matrixC2ESt5arrayImLm2EE
	.section	.text.unlikely,"ax",@progbits
.LCOLDB3:
	.section	.text.startup,"ax",@progbits
.LHOTB3:
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB1955:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA1955
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
	movl	$20, %esi
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movl	$21, %edx
	subq	$120, %rsp
	.cfi_def_cfa_offset 176
	leaq	16(%rsp), %rbp
	movq	%rbp, %rdi
	movq	%fs:40, %rax
	movq	%rax, 104(%rsp)
	xorl	%eax, %eax
.LEHB2:
	call	_ZN6matrixC1ESt5arrayImLm2EE
.LEHE2:
	leaq	64(%rsp), %rdi
	movl	$20, %esi
	movl	$21, %edx
.LEHB3:
	call	_ZN6matrixC1ESt5arrayImLm2EE
.LEHE3:
	movq	%rbp, %rdi
.LEHB4:
	call	_Z27initializeDisplacementArrayR6matrix@PLT
.LEHE4:
	movq	16(%rsp), %r15
	testq	%r15, %r15
	je	.L54
	movq	24(%rsp), %r12
	testq	%r12, %r12
	je	.L54
	vmovsd	.LC2(%rip), %xmm1
	vxorps	%xmm3, %xmm3, %xmm3
	leaq	0(,%r12,8), %rdx
	xorl	%r14d, %r14d
	xorl	%r13d, %r13d
	vxorpd	%xmm4, %xmm4, %xmm4
	.p2align 4
	.p2align 3
.L56:
	testq	%r13, %r13
	js	.L61
	vcvtsi2sdq	%r13, %xmm3, %xmm2
.L62:
	movq	32(%rsp), %rbp
	vmulsd	%xmm2, %xmm2, %xmm2
	xorl	%ebx, %ebx
	addq	%r14, %rbp
	.p2align 4
	.p2align 3
.L63:
	vcvtsi2sdq	%rbx, %xmm3, %xmm0
	vfmadd132sd	%xmm0, %xmm2, %xmm0
	vucomisd	%xmm0, %xmm4
	ja	.L93
	vsqrtsd	%xmm0, %xmm0, %xmm0
	vdivsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, 0(%rbp,%rbx,8)
	incq	%rbx
	cmpq	%r12, %rbx
	jne	.L63
.L60:
	incq	%r13
	addq	%rdx, %r14
	cmpq	%r15, %r13
	jne	.L56
.L54:
	movq	80(%rsp), %rdi
	testq	%rdi, %rdi
	je	.L64
	movq	96(%rsp), %rsi
	subq	%rdi, %rsi
	call	_ZdlPvm@PLT
.L64:
	movq	32(%rsp), %rdi
	testq	%rdi, %rdi
	je	.L65
	movq	48(%rsp), %rsi
	subq	%rdi, %rsi
	call	_ZdlPvm@PLT
.L65:
	movq	104(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L97
	addq	$120, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	xorl	%eax, %eax
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
.L61:
	.cfi_restore_state
	movq	%r13, %rax
	movq	%r13, %rcx
	shrq	%rax
	andl	$1, %ecx
	orq	%rcx, %rax
	vcvtsi2sdq	%rax, %xmm3, %xmm2
	vaddsd	%xmm2, %xmm2, %xmm2
	jmp	.L62
.L93:
	movq	%rdx, 8(%rsp)
	vmovsd	%xmm2, (%rsp)
	call	sqrt@PLT
	movq	.LC2(%rip), %rax
	vxorpd	%xmm4, %xmm4, %xmm4
	vxorps	%xmm3, %xmm3, %xmm3
	vmovsd	(%rsp), %xmm2
	movq	8(%rsp), %rdx
	vmovq	%rax, %xmm1
	vdivsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, 0(%rbp,%rbx,8)
	incq	%rbx
	cmpq	%r12, %rbx
	jne	.L63
	jmp	.L60
.L97:
	call	__stack_chk_fail@PLT
.L72:
	endbr64
	movq	%rax, %rbp
	jmp	.L66
.L71:
	endbr64
	movq	%rax, %rbp
	vzeroupper
	jmp	.L68
	.section	.gcc_except_table,"a",@progbits
.LLSDA1955:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE1955-.LLSDACSB1955
.LLSDACSB1955:
	.uleb128 .LEHB2-.LFB1955
	.uleb128 .LEHE2-.LEHB2
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB3-.LFB1955
	.uleb128 .LEHE3-.LEHB3
	.uleb128 .L71-.LFB1955
	.uleb128 0
	.uleb128 .LEHB4-.LFB1955
	.uleb128 .LEHE4-.LEHB4
	.uleb128 .L72-.LFB1955
	.uleb128 0
.LLSDACSE1955:
	.section	.text.startup
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDAC1955
	.type	main.cold, @function
main.cold:
.LFSB1955:
.L66:
	.cfi_def_cfa_offset 176
	.cfi_offset 3, -56
	.cfi_offset 6, -48
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	.cfi_offset 14, -24
	.cfi_offset 15, -16
	movq	80(%rsp), %rdi
	movq	96(%rsp), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	je	.L95
	vzeroupper
	call	_ZdlPvm@PLT
.L68:
	movq	32(%rsp), %rdi
	movq	48(%rsp), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	je	.L69
	call	_ZdlPvm@PLT
.L69:
	movq	%rbp, %rdi
.LEHB5:
	call	_Unwind_Resume@PLT
.LEHE5:
.L95:
	vzeroupper
	jmp	.L68
	.cfi_endproc
.LFE1955:
	.section	.gcc_except_table
.LLSDAC1955:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSEC1955-.LLSDACSBC1955
.LLSDACSBC1955:
	.uleb128 .LEHB5-.LCOLDB3
	.uleb128 .LEHE5-.LEHB5
	.uleb128 0
	.uleb128 0
.LLSDACSEC1955:
	.section	.text.unlikely
	.section	.text.startup
	.size	main, .-main
	.section	.text.unlikely
	.size	main.cold, .-main.cold
.LCOLDE3:
	.section	.text.startup
.LHOTE3:
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC2:
	.long	1841940611
	.long	1070882608
	.hidden	DW.ref.__gxx_personality_v0
	.weak	DW.ref.__gxx_personality_v0
	.section	.data.rel.local.DW.ref.__gxx_personality_v0,"awG",@progbits,DW.ref.__gxx_personality_v0,comdat
	.align 8
	.type	DW.ref.__gxx_personality_v0, @object
	.size	DW.ref.__gxx_personality_v0, 8
DW.ref.__gxx_personality_v0:
	.quad	__gxx_personality_v0
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
