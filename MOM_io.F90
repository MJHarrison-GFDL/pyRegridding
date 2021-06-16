!> Routines for error handling and I/O management
module MOM_io

  use MOM_string_functions, only : lowercase, slasher
  use MOM_error_handler, only : MOM_error, NOTE, FATAL, WARNING
implicit none ; private

public :: stdout,stderr
public :: verify_variable_units
public :: lowercase, slasher, file_exists, field_exists, field_size
public :: MOM_read_data, ensembler

integer :: stdout = 5
integer :: stderr = 6

contains

  !> Verify that a file contains a named variable with the expected units.
subroutine verify_variable_units(filename, varname, expected_units, msg, ierr, alt_units)
  character(len=*),           intent(in)    :: filename  !< File name
  character(len=*),           intent(in)    :: varname   !< Variable name
  character(len=*),           intent(in)    :: expected_units !< Expected units of variable
  character(len=*),           intent(inout) :: msg       !< Message to use for errors
  logical,                    intent(out)   :: ierr      !< True if an error occurs
  character(len=*), optional, intent(in)    :: alt_units !< Alterate acceptable units of variable

  ! Local variables
  character (len=200) :: units
  logical :: units_correct, success
  integer :: i, ncid, status, vid

  ierr = .false. ; msg = '' ; return


end subroutine verify_variable_units

logical function file_exists(filename)
  character(len=*),                 intent(in) :: filename   !< The name of the file being inquired about

  file_exists=.true.
end function file_exists

logical function field_exists(filename, fieldname)
  character(len=*),                 intent(in) :: filename   !< The name of the file being inquired about
  character(len=*),                 intent(in) :: fieldname !< The name of the field being sought

  field_exists=.true.
end function field_exists


subroutine field_size(filename, fieldname, siz)
  character(len=*),                 intent(in) :: filename   !< The name of the file being inquired about
  character(len=*),                 intent(in) :: fieldname !< The name of the field being sought
  integer, dimension(4) :: siz

  siz=-1
end subroutine field_size


  !> This routine uses the fms_io subroutine read_data to read a 1-D data field named
!! "fieldname" from a single or domain-decomposed file "filename".
subroutine MOM_read_data(filename, fieldname, data, timelevel, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before they are returned.


end subroutine MOM_read_data

!> Returns a name with "%#E" or "%E" replaced with the ensemble member number.
function ensembler(name, ens_no_in) result(en_nm)
  character(len=*),  intent(in) :: name       !< The name to be modified
  integer, optional, intent(in) :: ens_no_in  !< The number of the current ensemble member
  character(len=len(name)) :: en_nm  !< The name encoded with the ensemble number

  ! This function replaces "%#E" or "%E" with the ensemble number anywhere it
  ! occurs in name, with %E using 4 or 6 digits (depending on the ensemble size)
  ! and %#E using # digits, where # is a number from 1 to 9.

  character(len=len(name)) :: tmp
  character(10) :: ens_num_char
  character(3)  :: code_str
  integer :: ens_no
  integer :: n, is, ie

  en_nm = trim(name)
  if (index(name,"%") == 0) return

  if (present(ens_no_in)) then
    ens_no = ens_no_in
  else
    ens_no = 0
  endif

  write(ens_num_char, '(I10)') ens_no ; ens_num_char = adjustl(ens_num_char)
  do
    is = index(en_nm,"%E")
    if (is == 0) exit
    if (len(en_nm) < len(trim(en_nm)) + len(trim(ens_num_char)) - 2) &
      call MOM_error(FATAL, "MOM_io ensembler: name "//trim(name)// &
      " is not long enough for %E expansion for ens_no "//trim(ens_num_char))
    tmp = en_nm(1:is-1)//trim(ens_num_char)//trim(en_nm(is+2:))
    en_nm = tmp
  enddo

  if (index(name,"%") == 0) return

  write(ens_num_char, '(I10.10)') ens_no
  do n=1,9 ; do
    write(code_str, '("%",I1,"E")') n

    is = index(en_nm,code_str)
    if (is == 0) exit
    if (ens_no < 10**n) then
      if (len(en_nm) < len(trim(en_nm)) + n-3) call MOM_error(FATAL, &
        "MOM_io ensembler: name "//trim(name)//" is not long enough for %E expansion.")
      tmp = en_nm(1:is-1)//trim(ens_num_char(11-n:10))//trim(en_nm(is+3:))
    else
      call MOM_error(FATAL, "MOM_io ensembler: Ensemble number is too large "//&
          "to be encoded with "//code_str//" in "//trim(name))
    endif
    en_nm = tmp
  enddo ; enddo

end function ensembler

end module MOM_io
