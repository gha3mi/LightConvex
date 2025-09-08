program check
   ! Fortran Standard.
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   ! Unit-test utilities.
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   ! Collection of test problems.
   use TestLinalg, only: collect_linalg_tests
   use TestLinearPrograms, only: collect_dense_standard_simplex_problems
   implicit none(external)

! Unit-test related.
   integer :: status, i
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   ! Collection of test suites.
   status = 0
   testsuites = [new_testsuite("Linear Algebra", collect_linalg_tests)]
   testsuites = [testsuites, new_testsuite("Dense Standard Simplex", collect_dense_standard_simplex_problems)]

   ! Run all the test suites.
   do i = 1, size(testsuites)
      write (*, *) "-------------------------------"
      write (error_unit, fmt) "Testing :", testsuites(i)%name
      write (*, *) "-------------------------------"
      write (*, *)
      call run_testsuite(testsuites(i)%collect, error_unit, status)
   end do

   if (status > 0) then
      write (error_unit, '(i0, 1x, a)') status, "test(s) failed!"
   else if (status == 0) then
      write (output_unit, *) "All tests succesfully passed!", new_line('a')
   end if
end program check
