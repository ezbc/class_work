program hello_world_private

  use omp_lib
  implicit none

  integer :: id

  !$omp parallel num_threads(4) private(id)
  id = omp_get_thread_num()
  print *, 'Hello from thread', id
  !$omp end parallel

end program hello_world_private
