	.file	"MlmsWriteToFile.cpp"
	.text
	.section	.text._ZNKSt5ctypeIcE8do_widenEc,"axG",@progbits,_ZNKSt5ctypeIcE8do_widenEc,comdat
	.align 2
	.p2align 4
	.weak	_ZNKSt5ctypeIcE8do_widenEc
	.type	_ZNKSt5ctypeIcE8do_widenEc, @function
_ZNKSt5ctypeIcE8do_widenEc:
.LFB2185:
	.cfi_startproc
	endbr64
	movl	%esi, %eax
	ret
	.cfi_endproc
.LFE2185:
	.size	_ZNKSt5ctypeIcE8do_widenEc, .-_ZNKSt5ctypeIcE8do_widenEc
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC0:
	.string	"\t"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB2:
	.text
.LHOTB2:
	.p2align 4
	.globl	_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
	.type	_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE, @function
_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
.LFB3635:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA3635
	endbr64
	leaq	8(%rsp), %r10
	.cfi_def_cfa 10, 0
	andq	$-32, %rsp
	pushq	-8(%r10)
	pushq	%rbp
	movq	%rsp, %rbp
	.cfi_escape 0x10,0x6,0x2,0x76,0
	pushq	%r15
	pushq	%r14
	pushq	%r13
	.cfi_escape 0x10,0xf,0x2,0x76,0x78
	.cfi_escape 0x10,0xe,0x2,0x76,0x70
	.cfi_escape 0x10,0xd,0x2,0x76,0x68
	leaq	-344(%rbp), %r15
	pushq	%r12
	pushq	%r10
	.cfi_escape 0xf,0x3,0x76,0x58,0x6
	.cfi_escape 0x10,0xc,0x2,0x76,0x60
	pushq	%rbx
	subq	$576, %rsp
	.cfi_escape 0x10,0x3,0x2,0x76,0x50
	movq	%rdi, %r12
	movq	%r15, %rdi
	movq	%rsi, %rbx
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	movq	%r15, -608(%rbp)
	call	_ZNSt8ios_baseC2Ev@PLT
	leaq	16+_ZTVSt9basic_iosIcSt11char_traitsIcEE(%rip), %rax
	vpxor	%xmm0, %xmm0, %xmm0
	movq	%rax, -344(%rbp)
	movq	8+_ZTTSt14basic_ofstreamIcSt11char_traitsIcEE(%rip), %rax
	movw	$0, -120(%rbp)
	vmovdqa	%ymm0, -112(%rbp)
	leaq	-592(%rbp), %r13
	xorl	%esi, %esi
	movq	$0, -128(%rbp)
	movq	-24(%rax), %rdi
	movq	%rax, -592(%rbp)
	movq	16+_ZTTSt14basic_ofstreamIcSt11char_traitsIcEE(%rip), %rax
	addq	%r13, %rdi
	movq	%rax, (%rdi)
	vzeroupper
.LEHB0:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E@PLT
.LEHE0:
	leaq	24+_ZTVSt14basic_ofstreamIcSt11char_traitsIcEE(%rip), %rax
	movq	%rax, -592(%rbp)
	addq	$40, %rax
	movq	%rax, -344(%rbp)
	leaq	-584(%rbp), %rax
	movq	%rax, %rdi
	movq	%rax, %r14
	movq	%rax, -600(%rbp)
.LEHB1:
	call	_ZNSt13basic_filebufIcSt11char_traitsIcEEC1Ev@PLT
.LEHE1:
	movq	%r14, %rsi
	movq	%r15, %rdi
.LEHB2:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E@PLT
.LEHE2:
	movq	(%rbx), %rsi
	movq	-600(%rbp), %rdi
	movl	$48, %edx
.LEHB3:
	call	_ZNSt13basic_filebufIcSt11char_traitsIcEE4openEPKcSt13_Ios_Openmode@PLT
	movq	-592(%rbp), %rdx
	movq	-24(%rdx), %rdi
	addq	%r13, %rdi
	testq	%rax, %rax
	je	.L37
	xorl	%esi, %esi
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate@PLT
.L11:
	xorl	%r14d, %r14d
	cmpq	$0, (%r12)
	leaq	.LC0(%rip), %r15
	je	.L10
	.p2align 4
	.p2align 3
.L9:
	movq	8(%r12), %rax
	xorl	%ebx, %ebx
	testq	%rax, %rax
	je	.L16
	.p2align 4
	.p2align 3
.L12:
	movq	16(%r12), %rdx
	imulq	%r14, %rax
	movq	%r13, %rdi
	addq	%rbx, %rax
	vmovsd	(%rdx,%rax,8), %xmm0
	call	_ZNSo9_M_insertIdEERSoT_@PLT
	movl	$1, %edx
	movq	%r15, %rsi
	movq	%r13, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	8(%r12), %rax
	incq	%rbx
	cmpq	%rbx, %rax
	ja	.L12
.L16:
	movq	-592(%rbp), %rax
	movq	-24(%rax), %rax
	movq	-352(%rbp,%rax), %rbx
	testq	%rbx, %rbx
	je	.L38
	cmpb	$0, 56(%rbx)
	je	.L14
	movsbl	67(%rbx), %esi
.L15:
	movq	%r13, %rdi
	call	_ZNSo3putEc@PLT
	movq	%rax, %rdi
	call	_ZNSo5flushEv@PLT
	incq	%r14
	cmpq	%r14, (%r12)
	ja	.L9
.L10:
	movq	-600(%rbp), %rdi
	call	_ZNSt13basic_filebufIcSt11char_traitsIcEE5closeEv@PLT
.LEHE3:
	vmovq	.LC1(%rip), %xmm2
	leaq	16+_ZTVSt13basic_filebufIcSt11char_traitsIcEE(%rip), %rdx
	vpinsrq	$1, %rdx, %xmm2, %xmm1
	vmovdqa	%xmm1, -624(%rbp)
	testq	%rax, %rax
	je	.L39
.L17:
	vmovdqa	-624(%rbp), %xmm3
	leaq	64+_ZTVSt14basic_ofstreamIcSt11char_traitsIcEE(%rip), %rax
	movq	-600(%rbp), %rdi
	movq	%rax, -344(%rbp)
	vmovdqa	%xmm3, -592(%rbp)
.LEHB4:
	call	_ZNSt13basic_filebufIcSt11char_traitsIcEE5closeEv@PLT
.LEHE4:
.L19:
	leaq	-480(%rbp), %rdi
	call	_ZNSt12__basic_fileIcED1Ev@PLT
	leaq	16+_ZTVSt15basic_streambufIcSt11char_traitsIcEE(%rip), %rax
	leaq	-528(%rbp), %rdi
	movq	%rax, -584(%rbp)
	call	_ZNSt6localeD1Ev@PLT
	movq	8+_ZTTSt14basic_ofstreamIcSt11char_traitsIcEE(%rip), %rax
	movq	16+_ZTTSt14basic_ofstreamIcSt11char_traitsIcEE(%rip), %rcx
	movq	-608(%rbp), %rdi
	movq	%rax, -592(%rbp)
	movq	-24(%rax), %rax
	movq	%rcx, -592(%rbp,%rax)
	leaq	16+_ZTVSt9basic_iosIcSt11char_traitsIcEE(%rip), %rax
	movq	%rax, -344(%rbp)
	call	_ZNSt8ios_baseD2Ev@PLT
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L40
	addq	$576, %rsp
	popq	%rbx
	popq	%r10
	.cfi_remember_state
	.cfi_def_cfa 10, 0
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	leaq	-8(%r10), %rsp
	.cfi_def_cfa 7, 8
	ret
	.p2align 4
	.p2align 3
.L14:
	.cfi_restore_state
	movq	%rbx, %rdi
.LEHB5:
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%rbx), %rax
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rcx
	movl	$10, %esi
	movq	48(%rax), %rax
	cmpq	%rcx, %rax
	je	.L15
	movq	%rbx, %rdi
	call	*%rax
	movsbl	%al, %esi
	jmp	.L15
.L37:
	movl	32(%rdi), %esi
	orl	$4, %esi
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate@PLT
	jmp	.L11
.L39:
	movq	-592(%rbp), %rax
	movq	-24(%rax), %rdi
	addq	%r13, %rdi
	movl	32(%rdi), %esi
	orl	$4, %esi
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate@PLT
	jmp	.L17
.L40:
	call	__stack_chk_fail@PLT
.L38:
	call	_ZSt16__throw_bad_castv@PLT
.LEHE5:
.L26:
	endbr64
	movq	%rax, %r12
	vzeroupper
	jmp	.L7
.L29:
	endbr64
	movq	%rax, %rdi
	jmp	.L18
.L28:
	endbr64
	movq	%rax, %rbx
	jmp	.L5
.L25:
	endbr64
	movq	%rax, %r12
	jmp	.L20
.L27:
	endbr64
	movq	%rax, %rbx
	vzeroupper
	jmp	.L6
	.globl	__gxx_personality_v0
	.section	.gcc_except_table,"a",@progbits
	.align 4
.LLSDA3635:
	.byte	0xff
	.byte	0x9b
	.uleb128 .LLSDATT3635-.LLSDATTD3635
.LLSDATTD3635:
	.byte	0x1
	.uleb128 .LLSDACSE3635-.LLSDACSB3635
.LLSDACSB3635:
	.uleb128 .LEHB0-.LFB3635
	.uleb128 .LEHE0-.LEHB0
	.uleb128 .L26-.LFB3635
	.uleb128 0
	.uleb128 .LEHB1-.LFB3635
	.uleb128 .LEHE1-.LEHB1
	.uleb128 .L27-.LFB3635
	.uleb128 0
	.uleb128 .LEHB2-.LFB3635
	.uleb128 .LEHE2-.LEHB2
	.uleb128 .L28-.LFB3635
	.uleb128 0
	.uleb128 .LEHB3-.LFB3635
	.uleb128 .LEHE3-.LEHB3
	.uleb128 .L25-.LFB3635
	.uleb128 0
	.uleb128 .LEHB4-.LFB3635
	.uleb128 .LEHE4-.LEHB4
	.uleb128 .L29-.LFB3635
	.uleb128 0x1
	.uleb128 .LEHB5-.LFB3635
	.uleb128 .LEHE5-.LEHB5
	.uleb128 .L25-.LFB3635
	.uleb128 0
.LLSDACSE3635:
	.byte	0x1
	.byte	0
	.align 4
	.long	0

.LLSDATT3635:
	.text
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDAC3635
	.type	_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE.cold, @function
_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE.cold:
.LFSB3635:
.L5:
	.cfi_escape 0xf,0x3,0x76,0x58,0x6
	.cfi_escape 0x10,0x3,0x2,0x76,0x50
	.cfi_escape 0x10,0x6,0x2,0x76,0
	.cfi_escape 0x10,0xc,0x2,0x76,0x60
	.cfi_escape 0x10,0xd,0x2,0x76,0x68
	.cfi_escape 0x10,0xe,0x2,0x76,0x70
	.cfi_escape 0x10,0xf,0x2,0x76,0x78
	movq	-600(%rbp), %rdi
	vzeroupper
	call	_ZNSt13basic_filebufIcSt11char_traitsIcEED1Ev@PLT
.L6:
	movq	8+_ZTTSt14basic_ofstreamIcSt11char_traitsIcEE(%rip), %rax
	movq	16+_ZTTSt14basic_ofstreamIcSt11char_traitsIcEE(%rip), %rcx
	movq	%rbx, %r12
	movq	%rax, -592(%rbp)
	movq	-24(%rax), %rax
	movq	%rcx, -592(%rbp,%rax)
.L7:
	movq	-608(%rbp), %rdi
	leaq	16+_ZTVSt9basic_iosIcSt11char_traitsIcEE(%rip), %rax
	movq	%rax, -344(%rbp)
	call	_ZNSt8ios_baseD2Ev@PLT
	movq	%r12, %rdi
.LEHB6:
	call	_Unwind_Resume@PLT
.LEHE6:
.L18:
	vzeroupper
	call	__cxa_begin_catch@PLT
	call	__cxa_end_catch@PLT
	jmp	.L19
.L20:
	movq	%r13, %rdi
	vzeroupper
	call	_ZNSt14basic_ofstreamIcSt11char_traitsIcEED1Ev@PLT
	movq	%r12, %rdi
.LEHB7:
	call	_Unwind_Resume@PLT
.LEHE7:
	.cfi_endproc
.LFE3635:
	.section	.gcc_except_table
	.align 4
.LLSDAC3635:
	.byte	0xff
	.byte	0x9b
	.uleb128 .LLSDATTC3635-.LLSDATTDC3635
.LLSDATTDC3635:
	.byte	0x1
	.uleb128 .LLSDACSEC3635-.LLSDACSBC3635
.LLSDACSBC3635:
	.uleb128 .LEHB6-.LCOLDB2
	.uleb128 .LEHE6-.LEHB6
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB7-.LCOLDB2
	.uleb128 .LEHE7-.LEHB7
	.uleb128 0
	.uleb128 0
.LLSDACSEC3635:
	.byte	0x1
	.byte	0
	.align 4
	.long	0

.LLSDATTC3635:
	.section	.text.unlikely
	.text
	.size	_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE, .-_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
	.section	.text.unlikely
	.size	_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE.cold, .-_Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE.cold
.LCOLDE2:
	.text
.LHOTE2:
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.type	_GLOBAL__sub_I__Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE, @function
_GLOBAL__sub_I__Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
.LFB4461:
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
.LFE4461:
	.size	_GLOBAL__sub_I__Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE, .-_GLOBAL__sub_I__Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I__Z11writeToFileRK6matrixNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.section	.data.rel.ro,"aw"
	.align 8
.LC1:
	.quad	_ZTVSt14basic_ofstreamIcSt11char_traitsIcEE+24
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
