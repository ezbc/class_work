program square_nums

  use omp_lib
  implicit none

  integer :: i
  integer :: j(4)

  !$omp parallel num_threads(4)
  i = omp_get_thread_num() + 1
  j(i) = i**2
  !$omp end parallel

  print *,'Square numbers:', j

end program square_nums
