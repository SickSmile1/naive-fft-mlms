	.file	"MlmsTimer.cpp"
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
.LFB4027:
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
.LFE4027:
	.size	_ZNSt6vectorIdSaIdEE17_M_default_appendEm, .-_ZNSt6vectorIdSaIdEE17_M_default_appendEm
	.section	.text.unlikely,"ax",@progbits
.LCOLDB8:
	.text
.LHOTB8:
	.p2align 4
	.globl	_Z13runTimerLoopsv
	.type	_Z13runTimerLoopsv, @function
_Z13runTimerLoopsv:
.LFB3550:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA3550
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
	movl	$40, %esi
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$216, %rsp
	.cfi_def_cfa_offset 272
	vmovdqa	.LC2(%rip), %xmm0
	movq	%fs:40, %rax
	movq	%rax, 200(%rsp)
	xorl	%eax, %eax
	leaq	80(%rsp), %rdi
	movq	$0, 96(%rsp)
	vmovdqa	%xmm0, 64(%rsp)
	vpxor	%xmm0, %xmm0, %xmm0
	vmovdqa	%xmm0, 80(%rsp)
.LEHB0:
	call	_ZNSt6vectorIdSaIdEE17_M_default_appendEm
.LEHE0:
	xorl	%r12d, %r12d
	movl	$480, %ebx
	leaq	128(%rsp), %r15
.L41:
	vmovq	%rbx, %xmm2
	movq	%rbx, %rbp
	movq	%r15, %rdi
	movq	$0, 144(%rsp)
	vpunpcklqdq	%xmm2, %xmm2, %xmm7
	vxorpd	%xmm2, %xmm2, %xmm2
	imulq	%rbx, %rbp
	vcvtsi2sdq	%rbx, %xmm2, %xmm0
	vmovsd	.LC3(%rip), %xmm2
	vmovsd	%xmm0, (%rsp)
	vmovdqa	%xmm7, 16(%rsp)
	vmovdqa	%xmm7, 112(%rsp)
	movq	%rbp, %rsi
	vdivsd	%xmm0, %xmm2, %xmm3
	vpxor	%xmm0, %xmm0, %xmm0
	vmovdqa	%xmm0, 128(%rsp)
	vmovsd	%xmm3, 8(%rsp)
.LEHB1:
	call	_ZNSt6vectorIdSaIdEE17_M_default_appendEm
.LEHE1:
	vmovdqa	16(%rsp), %xmm7
	leaq	176(%rsp), %rax
	vpxor	%xmm0, %xmm0, %xmm0
	movq	%rbp, %rsi
	movq	%rax, %rdi
	vmovdqa	%xmm0, 176(%rsp)
	movq	$0, 192(%rsp)
	movq	%rax, 16(%rsp)
	vmovdqa	%xmm7, 160(%rsp)
.LEHB2:
	call	_ZNSt6vectorIdSaIdEE17_M_default_appendEm
.LEHE2:
	vmovsd	8(%rsp), %xmm6
	leaq	160(%rsp), %r13
	vmovsd	.LC4(%rip), %xmm3
	vmovsd	.LC3(%rip), %xmm7
	movq	%r13, %rdi
	movq	.LC6(%rip), %rax
	vmovq	%rax, %xmm2
	vdivsd	%xmm6, %xmm3, %xmm1
	vmovsd	.LC5(%rip), %xmm3
	vmulsd	.LC5(%rip), %xmm1, %xmm1
	vdivsd	%xmm6, %xmm7, %xmm0
	vfmsub132sd	%xmm0, %xmm1, %xmm3
	vfmadd231sd	.LC5(%rip), %xmm0, %xmm1
	vmovsd	%xmm3, %xmm3, %xmm0
.LEHB3:
	call	_Z23initializePressureArrayR6matrixddd@PLT
	leaq	112(%rsp), %r14
	movq	%r14, %rdi
	call	_Z27initializeDisplacementArrayR6matrix@PLT
	call	_ZNSt6chrono3_V212steady_clock3nowEv@PLT
	vmovsd	8(%rsp), %xmm0
	vxorpd	%xmm1, %xmm1, %xmm1
	movq	%r13, %rsi
	movq	%rax, %rbp
	movq	.LC6(%rip), %rax
	movq	%r14, %rdi
	vmovq	%rax, %xmm2
	call	_Z16calculation_loopR6matrixRKS_ddd@PLT
.LEHE3:
	call	_ZNSt6chrono3_V212steady_clock3nowEv@PLT
	vxorpd	%xmm4, %xmm4, %xmm4
	vmovsd	(%rsp), %xmm5
	addq	$20, %rbx
	subq	%rbp, %rax
	movq	72(%rsp), %rdx
	movq	176(%rsp), %rdi
	vcvtsi2sdq	%rax, %xmm4, %xmm0
	vdivsd	.LC7(%rip), %xmm0, %xmm0
	movq	80(%rsp), %rax
	imulq	%r12, %rdx
	vunpcklpd	%xmm0, %xmm5, %xmm0
	vmovupd	%xmm0, (%rax,%rdx,8)
	testq	%rdi, %rdi
	je	.L52
	movq	192(%rsp), %rsi
	subq	%rdi, %rsi
	call	_ZdlPvm@PLT
	.p2align 4
	.p2align 3
.L52:
	movq	128(%rsp), %rdi
	testq	%rdi, %rdi
	je	.L53
	movq	144(%rsp), %rsi
	incq	%r12
	subq	%rdi, %rsi
	call	_ZdlPvm@PLT
	cmpq	$580, %rbx
	jne	.L41
.L55:
	movq	16(%rsp), %rbx
	movabsq	$3347132671979252084, %rax
	leaq	64(%rsp), %rdi
	movq	%r13, %rsi
	movq	%rax, 176(%rsp)
	movw	$30836, 184(%rsp)
	movb	$116, 186(%rsp)
	movq	$11, 168(%rsp)
	movb	$0, 187(%rsp)
	movq	%rbx, 160(%rsp)
.LEHB4:
	call	_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE4:
	movq	160(%rsp), %rdi
	cmpq	%rbx, %rdi
	je	.L56
	movq	176(%rsp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L56:
	movq	80(%rsp), %rdi
	testq	%rdi, %rdi
	je	.L40
	movq	96(%rsp), %rsi
	subq	%rdi, %rsi
	call	_ZdlPvm@PLT
.L40:
	movq	200(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L112
	addq	$216, %rsp
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
.L53:
	.cfi_restore_state
	incq	%r12
	cmpq	$580, %rbx
	jne	.L41
	jmp	.L55
.L112:
	call	__stack_chk_fail@PLT
.L67:
	endbr64
	movq	%rax, %rbp
	jmp	.L42
.L66:
	endbr64
	movq	%rax, %rbp
	jmp	.L61
.L69:
	endbr64
	movq	%rax, %rbp
	jmp	.L58
.L70:
	endbr64
	movq	%rax, %rbp
	jmp	.L58
.L68:
	endbr64
	movq	%rax, %rbp
	jmp	.L45
	.globl	__gxx_personality_v0
	.section	.gcc_except_table,"a",@progbits
.LLSDA3550:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE3550-.LLSDACSB3550
.LLSDACSB3550:
	.uleb128 .LEHB0-.LFB3550
	.uleb128 .LEHE0-.LEHB0
	.uleb128 .L67-.LFB3550
	.uleb128 0
	.uleb128 .LEHB1-.LFB3550
	.uleb128 .LEHE1-.LEHB1
	.uleb128 .L68-.LFB3550
	.uleb128 0
	.uleb128 .LEHB2-.LFB3550
	.uleb128 .LEHE2-.LEHB2
	.uleb128 .L69-.LFB3550
	.uleb128 0
	.uleb128 .LEHB3-.LFB3550
	.uleb128 .LEHE3-.LEHB3
	.uleb128 .L70-.LFB3550
	.uleb128 0
	.uleb128 .LEHB4-.LFB3550
	.uleb128 .LEHE4-.LEHB4
	.uleb128 .L66-.LFB3550
	.uleb128 0
.LLSDACSE3550:
	.text
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDAC3550
	.type	_Z13runTimerLoopsv.cold, @function
_Z13runTimerLoopsv.cold:
.LFSB3550:
.L61:
	.cfi_def_cfa_offset 272
	.cfi_offset 3, -56
	.cfi_offset 6, -48
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	.cfi_offset 14, -24
	.cfi_offset 15, -16
	movq	160(%rsp), %rdi
	cmpq	16(%rsp), %rdi
	je	.L105
	movq	176(%rsp), %rax
	leaq	1(%rax), %rsi
.L109:
	vzeroupper
.L110:
	call	_ZdlPvm@PLT
.L47:
	movq	80(%rsp), %rdi
	movq	96(%rsp), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	je	.L63
	call	_ZdlPvm@PLT
.L63:
	movq	%rbp, %rdi
.LEHB5:
	call	_Unwind_Resume@PLT
.L45:
	movq	128(%rsp), %rdi
	movq	144(%rsp), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	jne	.L109
.L105:
	vzeroupper
	jmp	.L47
.L42:
	movq	80(%rsp), %rdi
	movq	96(%rsp), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	je	.L101
	vzeroupper
	call	_ZdlPvm@PLT
.L43:
	movq	%rbp, %rdi
	call	_Unwind_Resume@PLT
.LEHE5:
.L58:
	movq	176(%rsp), %rdi
	movq	192(%rsp), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	je	.L104
	vzeroupper
	call	_ZdlPvm@PLT
.L51:
	movq	128(%rsp), %rdi
	movq	144(%rsp), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	jne	.L110
	jmp	.L47
.L101:
	vzeroupper
	jmp	.L43
.L104:
	vzeroupper
	jmp	.L51
	.cfi_endproc
.LFE3550:
	.section	.gcc_except_table
.LLSDAC3550:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSEC3550-.LLSDACSBC3550
.LLSDACSBC3550:
	.uleb128 .LEHB5-.LCOLDB8
	.uleb128 .LEHE5-.LEHB5
	.uleb128 0
	.uleb128 0
.LLSDACSEC3550:
	.section	.text.unlikely
	.text
	.size	_Z13runTimerLoopsv, .-_Z13runTimerLoopsv
	.section	.text.unlikely
	.size	_Z13runTimerLoopsv.cold, .-_Z13runTimerLoopsv.cold
.LCOLDE8:
	.text
.LHOTE8:
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.type	_GLOBAL__sub_I__Z13runTimerLoopsv, @function
_GLOBAL__sub_I__Z13runTimerLoopsv:
.LFB4256:
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
.LFE4256:
	.size	_GLOBAL__sub_I__Z13runTimerLoopsv, .-_GLOBAL__sub_I__Z13runTimerLoopsv
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I__Z13runTimerLoopsv
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC2:
	.quad	20
	.quad	2
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC3:
	.long	0
	.long	1083129856
	.align 8
.LC4:
	.long	0
	.long	1082081280
	.align 8
.LC5:
	.long	0
	.long	1071644672
	.align 8
.LC6:
	.long	0
	.long	1072693248
	.align 8
.LC7:
	.long	0
	.long	1104006501
	.hidden	DW.ref.__gxx_personality_v0
	.weak	DW.ref.__gxx_personality_v0
	.section	.data.rel.local.DW.ref.__gxx_personality_v0,"awG",@progbits,DW.ref.__gxx_personality_v0,comdat
	.align 8
	.type	DW.ref.__gxx_personality_v0, @object
	.size	DW.ref.__gxx_personality_v0, 8
DW.ref.__gxx_personality_v0:
	.quad	__gxx_personality_v0
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
