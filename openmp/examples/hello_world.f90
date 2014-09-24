program hello_world

  use omp_lib
  implicit none

  !$omp parallel num_threads(4)
  print *,'Hello from thread', omp_get_thread_num()
  !$omp end parallel

end program hello_world
