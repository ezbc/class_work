program hello_world
    use omp_lib
    implicit none
    integer :: id
    !$omp parallel private(id)
        id = omp_get_thread_num()
        print *, 'Hello from thread', id
    !$omp end parallel
end program hello_world

