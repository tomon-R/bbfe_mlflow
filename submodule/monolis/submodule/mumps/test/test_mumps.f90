program test_mumps
use, intrinsic:: iso_fortran_env, only: output_unit, stderr=>error_unit, int32

implicit none

include 'mpif.h'
include 'dmumps_struc.h'  ! per MUMPS manual
type(DMUMPS_STRUC) :: mumps_par

integer :: num_mpi
integer(int32) :: ierr

call mpi_init(ierr)
if (ierr /= 0) error stop 'MPI init error'

call MPI_COMM_size(MPI_COMM_WORLD, num_mpi, ierr)
if(ierr/=0) error stop 'problem getting number of MPI processes'
print '(a,i0,a)', 'using ',num_mpi,' MPI processes'

mumps_par%COMM = MPI_COMM_WORLD

mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call DMUMPS(mumps_par)

! must set ICNTL after initialization Job= -1 above

mumps_par%icntl(1) = stderr  ! error messages
mumps_par%icntl(2) = output_unit !  diagnostic, statistics, and warning messages
mumps_par%icntl(3) = output_unit ! global info, for the host (myid==0)
mumps_par%icntl(4) = 1           ! default is 2, this reduces verbosity

IF (mumps_par%INFOG(1) < 0) THEN
  WRITE(stderr,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
  "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
  "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)

  error stop
END IF

call mpi_finalize(ierr)
if (ierr /= 0) error stop 'MPI finalize error'

print *, 'OK: Mumps params'

end program
