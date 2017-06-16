!============<Finding const of the task>==================

subroutine find_const(g_a, g_b, g_c, g_m, m_c)
! приближенные значения углов ga,gb,gc,gm,go,  модуля  и аргумента С
use mod
integer(4), parameter :: n = 5                  ! число уравнений    
integer(4) :: itmax = 200                       ! кол-во итераций
real(8) g_a, g_b, g_c, g_m, m_c                 !
real(8) x(n), xguess(n), errrel, fnorm
complex(8) test
external nles
xguess = (/g_a, g_b, g_c, g_m, m_c/)            ! начальное приближение
errrel = eps                                    ! погрешность

call dneqnf(nles, errrel, n, itmax, xguess, x, fnorm)
call dwrrrn('x', 1, n, x, 1, 0)
print *, 'fnorm =', fnorm 

ga = x(1)
gb = x(2)
gc = x(3)
gm = x(4)
mod_C = x(5)

!! проверка метода ньютона, не совсем корректно
test = - cdexp(ii * gm) + (d2 - delta) * cdexp(ii * gc) + c1 + (delta - d1) * cdexp(ii * gb) + cdexp(ii * ga)
print *, 'test = ', test

arg_C = -d5 * (ga + (delta - d1) * gb + (d2 - delta)* gc - gm) + (3 / d2 - delta) * pi
u_inf = mod_C * v_inf
beta = alpha - arg_C
go = d2 * beta - gc + gm + pi
print *, alpha

Q = d2 * pi * u_inf * (dcos(gm - beta) - dcos(go - beta) - dcos(gc - beta) - dcos(beta))
Gamma = d2 * pi * u_inf * (dsin(gm - beta) - dsin(go - beta) - dsin(gc - beta) + dsin(beta))
! характерный размер b = l_1 + l_2
Cy = d2 * d2 * pi * mod_C * (dsin(gm - beta) - dsin(go - beta) - dsin(gc - beta) + dsin(beta)) / (l_1 + l_2)

end subroutine

subroutine nles(x, f, n)
! non-linear equations system
! подпрограмма вычисления функции
! следует обратить внимание на порядок параметров n, x, f было раньше
use mod
integer(4) n 
real(8) x(n), f(n), res3, res4, res5, ds_dgamma, errest3, errest4, errest5
external ds_dgamma
ga = x(1)
gb = x(2)
gc = x(3)
gm = x(4)
mod_C = x(5)

f(1) = dcos(x(1)) + (delta - d1) * dcos(x(2)) + d1 + (d2 - delta) * dcos(x(3)) - dcos(x(4))
f(2) = dsin(x(1)) + (delta - d1) * dsin(x(2)) + (d2 - delta) * dsin(x(3)) - dsin(x(4))

call dqdags(ds_dgamma, gb, ga, d0, eps, res3, errest)
f(3) = res3 - l_1

call dqdags(ds_dgamma, d0, gb, d0, eps, res4, errest)
f(4) = res4 - l_2

call dqdags(ds_dgamma, ga, gc, d0, eps, res5, errest)
f(5) = res5 - l_1 + l_3
end subroutine 

real(8) function ds_dgamma(g)
! g - gamma
use mod
real(8) g
ds_dgamma = d2*d2*mod_C * dabs(dsin((g - ga) / d2) * dsin(g / d2) / dsin((g - gm) / d2)) * dabs(dsin((g - gb) / d2))**(delta-d1) * dabs(dsin((g - gc) / d2))**(d2-delta) 
end function
