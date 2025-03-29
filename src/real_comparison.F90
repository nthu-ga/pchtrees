module real_comparison
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none

    ! APC: must admit that I used ChatGPT for this code.

    private  ! Hide everything except the public interface
    public :: real_equal

    interface real_equal
        module procedure real_equal_r4, real_equal_r8, real_equal_mixed
    end interface real_equal

contains

    ! Single-precision (real(4)) version
    function real_equal_r4(a, b, abs_tol, rel_tol) result(is_equal)
        real(real32), intent(in) :: a, b
        real(real32), intent(in), optional :: abs_tol, rel_tol
        logical :: is_equal
        real(real32) :: epsilon_abs, epsilon_rel, diff, scale

        ! Default tolerances based on machine precision
        epsilon_abs = merge(abs_tol, 1.0e-6_real32, present(abs_tol))
        epsilon_rel = merge(rel_tol, 1.0e-5_real32, present(rel_tol))

        diff = abs(a - b)
        scale = max(abs(a), abs(b), 1.0_real32)

        is_equal = (diff < epsilon_abs) .or. (diff < epsilon_rel * scale)
    end function real_equal_r4

    ! Double-precision (real(8)) version
    function real_equal_r8(a, b, abs_tol, rel_tol) result(is_equal)
        real(real64), intent(in) :: a, b
        real(real64), intent(in), optional :: abs_tol, rel_tol
        logical :: is_equal
        real(real64) :: epsilon_abs, epsilon_rel, diff, scale

        ! Default tolerances based on machine precision
        epsilon_abs = merge(abs_tol, 1.0e-12_real64, present(abs_tol))
        epsilon_rel = merge(rel_tol, 1.0e-10_real64, present(rel_tol))

        diff = abs(a - b)
        scale = max(abs(a), abs(b), 1.0_real64)

        is_equal = (diff < epsilon_abs) .or. (diff < epsilon_rel * scale)
    end function real_equal_r8

    ! Mixed-precision (real(4) and real(8)) version
    function real_equal_mixed(a, b, abs_tol, rel_tol) result(is_equal)
        real(real32), intent(in) :: a
        real(real64), intent(in) :: b
        real(real64), intent(in), optional :: abs_tol, rel_tol
        logical :: is_equal
        real(real64) :: a_dbl, epsilon_abs, epsilon_rel, diff, scale

        ! Promote single-precision input to double
        a_dbl = real(a, kind=real64)

        ! Default tolerances based on double precision
        epsilon_abs = merge(abs_tol, 1.0e-12_real64, present(abs_tol))
        epsilon_rel = merge(rel_tol, 1.0e-10_real64, present(rel_tol))

        diff = abs(a_dbl - b)
        scale = max(abs(a_dbl), abs(b), 1.0_real64)

        is_equal = (diff < epsilon_abs) .or. (diff < epsilon_rel * scale)
    end function real_equal_mixed

end module real_comparison
