program main
	print *
	!call test_constructor()
	!print *
	!call test_accessors()
	!print *
	call test_read_txt()
	print *
	call test_write_txt()
	print *
end program main

subroutine test_constructor()
	use, intrinsic :: iso_fortran_env, only: output_unit
	use element_mod
	implicit none

	character(len=2),dimension(:),allocatable :: symbol
	type(element) :: obj

	print *, "test_constructor()"
	print *, "obj = element((/'Fe','C '/),2)"
	print *, "obj%write_txt_formatted()"
	symbol = (/'Fe','C '/)
	obj = element(symbol,2)
	call obj%write_txt_formatted(unit=output_unit)
end subroutine test_constructor

subroutine test_accessors()
	!use element_mod
	!implicit none

	!type(element) :: obj

	!print *, "test_accessors()"

	!print *, "obj = element('C')"
	!print *, "obj%set_symbol('carbon')"
	!print *, "obj%print()"
	!obj = element('C')
	!call obj%set_symbol('carbon')
	!call obj%print()

	!print *, "obj = element('C')"
	!print *, "obj%get_symbol() ="
	!obj = element('C')
	!print *, obj%get_symbol()
end subroutine test_accessors

subroutine test_read_txt()
	use, intrinsic :: iso_fortran_env, only: output_unit
	use element_mod
	implicit none

	type(element) :: obj

	print *, "test_read_txt()"
	print *, "obj%read_txt()"
	print *, "obj%write_txt_formatted()"
	call obj%read_txt()
	call obj%write_txt_formatted(unit=output_unit)
end subroutine test_read_txt

subroutine test_write_txt()
	use element_mod
	implicit none

	character(len=2),dimension(:),allocatable :: symbol
	type(element) :: obj

	print *, "test_write_txt()"
	print *, "obj = element((/'Fe','C '/),2)"
	print *, "obj%write_txt()"
	symbol = (/'Fe','C '/)
	obj = element(symbol,2)
	call obj%write_txt()
end subroutine test_write_txt

