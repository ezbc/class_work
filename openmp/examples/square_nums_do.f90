program square_nums_do

  use omp_lib
  implicit none

  integer :: i, j(13)

  !$omp parallel num_threads(4)
  !$omp do
  do i = 1, SIZE(j)
    j(i) = i**2
  enddo
  !$omp end parallel

  print *,'Square numbers:', j

end program square_nums_do
