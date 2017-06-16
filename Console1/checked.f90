!����������� ��� ��� �� ����������, ������ ��� ������ �� ������, �� ������� � ������ 
    
!==========
subroutine test_newtone
use gen_mod
integer(4), parameter :: n = 4 ! ����� ���������    
integer(4) :: itmax = 200 ! ���-�� ��������
real(8) x(n), xguess(n), errrel, fnorm
external fcn 
xguess = d1 ! ��������� �����������
errrel = eps ! �����������
call dneqnf(fcn, errrel, n, itmax, xguess, x, fnorm)
call dwrrrn('x', 1, n, x, 1, 0)
print *, 'fnorm =', fnorm 
end subroutine 

subroutine fcn(x, f, n)
! ������������ ���������� �������
! ������� �������� �������� �� ������� ���������� n, x, f ���� ������
use gen_mod
integer(4) n 
real(8) x(n), f(n)
f(1) = x(1) + d2*x(2) + x(3) + d2*d2*x(4) - 20.7d0
f(2) = x(1) * x(1) + d2 * x(1) * x(2) + x(4) * x(4) * x(4) - 15.88d0
f(3) = x(1) * x(1) * x(1) + x(3) * x(3) + x(4) - 21.218d0
f(4) = 3.0d0 * x(2) + x(3) * x(4) - 21.1d0
end subroutine 

!==========
subroutine test_integ
!�������� ���������� ���������
use gen_mod
real(8) a, b, res, f, errest ! res - ��������� ��������������, errest - ������ ������. �������� ��������
external f
a = d0; b = d1 ! ������� ��������������
call dqdags(f, a, b, d0, eps, res, errest)
print *, res, errest
end subroutine 
real function f(x)
! ��� �������� ������ dqdags
use gen_mod
real(8) x
f = dlog(x) / dsqrt(x)
end
!==========