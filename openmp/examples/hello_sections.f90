program hello_sections

  use omp_lib
  implicit none

  !$omp parallel schedule (static)
  !$omp sections
  !$omp section
  print *,'Hello from thread', omp_get_thread_num()
  !$omp section
  print *,'Hola from thread ', omp_get_thread_num()
  !$omp section
  print *,'Aloha from thread', omp_get_thread_num()
  !$omp end sections
  !$omp end parallel

end program hello_sections
