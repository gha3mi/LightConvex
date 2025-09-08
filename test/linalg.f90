module TestLinalg
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_math, only: all_close
   use lightconvex_constants, only: ilp, dp
   use lightconvex_linalg, only: qr_type, qr, factmv
   implicit none(external)
   private

   public :: collect_linalg_tests
contains
   subroutine collect_linalg_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("QR testsuite", test_qr_testsuite)]
   end subroutine collect_linalg_tests

   subroutine test_qr_testsuite(error)
      type(error_type), allocatable, intent(out) :: error

      !----- QR matrix-vector product -----
      qr_matvec: block
         integer(ilp), parameter :: m = 10, n = 5
         real(dp) :: A(m, n), x(n), y(m), y_(m)
         type(qr_type) :: F

         !> Create random matrix and vector.
         call random_number(A); call random_number(x)

         !> Factorize matrix.
         F = qr(A)

         !> Reference matvec.
         y = matmul(A, x)

         !> Factorized matvec.
         y_ = 0.0_dp; call factmv(F, x, y_)

         call check(error, all_close(y, y_))
         if (allocated(error)) return
      end block qr_matvec

      transposed_qr_matvec: block
         integer(ilp), parameter :: m = 10, n = 5
         real(dp) :: A(m, n), x(m), y(n), y_(n)
         type(qr_type) :: F

         !> Create random matrix.
         call random_number(A); call random_number(x)

         !> Factorize matrix.
         F = qr(A)

         !> Reference matvec.
         y = matmul(transpose(A), x)

         !> Factorized matvec.
         y_ = 0.0_dp; call factmv(F, x, y_, trans="T")

         call check(error, all_close(y, y_))
         if (allocated(error)) return
      end block transposed_qr_matvec
   end subroutine test_qr_testsuite
end module TestLinalg
