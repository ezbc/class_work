program square_nums_pardo

  use omp_lib
  implicit none

  integer :: i, j(13)

  !$omp parallel do num_threads(4)
  do i = 1, SIZE(j)
    j(i) = i**2
  enddo
  !$omp end parallel do

  print *,'Square numbers:', j

end program square_nums_pardo
