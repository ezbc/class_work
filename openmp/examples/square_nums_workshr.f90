program square_nums_workshr

  use omp_lib
  implicit none

  integer :: i, j(13)

  j = [(i,i=1,SIZE(j))]

  !$omp parallel workshare num_threads(4)
  j = j**2
  !$omp end parallel workshare

  print *,'Square numbers:', j

end program square_nums_workshr
