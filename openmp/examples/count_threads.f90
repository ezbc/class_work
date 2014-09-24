program count_threads

  use omp_lib
  implicit none

  integer :: n

  n = 0

  !$omp parallel num_threads(64)
  n = n + 1
  !$omp end parallel

  print *,'There were', n, 'threads'

end program count_threads
