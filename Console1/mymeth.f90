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

arg_C = - d5 * (ga + (delta - d1) * gb + (d2 - delta) * gc - gm) + (d2 - delta) * pi 
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

function dz_real(g)
! реальная часть подынтегральной ф-ции
! g - гамма угол между 
use mod
real(8) g, dz_real, g_a, g_b, g_c, g_m, g_0
complex(8) zzeta
g_a = ga
g_b = gb
g_c = gc
g_m = gm
g_0 = d0
!dz_real = dreal(с * zzeta(g, ga) * zzeta(g, gb) ** (delta - d1) * zzeta(g, d0) * zzeta(g, gc) ** (2 - delta) / cdexp(ii * g * d2) / zzeta(g, gm))
dz_real = dreal( mod_C * cdexp(ii * arg_C) * zzeta(g, g_a) * (zzeta(g, g_b) ** (delta - d1)) * zzeta(g, g_0) * (zzeta(g, g_c) ** (2 - delta)) / cdexp(ii * g * d2) / zzeta(g, g_m))
end function

function zzeta(x, y)
! ф-ция написана для сокращения кода при вычислении dz(g)
! разложение (zeta - zeta0) - с учетом выбора ветвей
use mod
real(8) x, y
complex(8) zzeta
zzeta = cdexp(ii * (x + y + pi * dsign(d1, x - y)) / d2) * dabs(d2 * dsin((x - y) / d2))
!zzeta = cdexp(ii * (x + y + pi) / d2) * dabs(d2 * dsin((x - y) / d2))
end function 

subroutine write_complex_array_to_file(array, n, file_name, str_length, port)
! запись в файл массива компексных чисел
! array - массив компексных чисел
! n - кол-во элементов массива
! file_name - имя файла
! str_length - кол-во символов 
! port - номер порта вывода
integer(4) i, n, str_length, port 
complex(8) array(n)
character(str_length) file_name
open(port, FILE = file_name) 
    write(port,*) ' VARIABLES = "X", "Y" '
    do i = 1, n
        write(port,"(F12.5,' ', F12.5)") dreal(array(i)), dimag(array(i))
    end do
close(port)
port = port + 1
end subroutine

function sgu(x)
! сгущает точки к середине
use mod
real(8) x, sgu
sgu = (x - d5) * (x - d5) * (x - d5) / (d5 * d5) + d5
!sgu=d1 - ( atan( tan(( (1-x)*(1-x) )*pi/d2) ) )*d2/pi
end function    

function dz_dzeta(x)
! поынтегральная ф.
! упрощал эту функцию, поэтому возможны ошибки 09.04.2017
use mod
real(4) x
complex(8) C, dz_dzeta
dz_dzeta = mod_C * dabs(d2*dsin(d5*(x-ga))) * (dabs(d2*dsin(d5*(x-gb))) ** (delta - d1)) * dabs(d2*dsin(d5*x)) * (dabs(d2*dsin(d5*(x-gc))) ** (d2 - delta)) / dabs(d2*dsin(d5*(x-gm))) * cdexp(ii*(arg_C - g + d5*(ga + (delta - d1)*gb + (d2 - delta)*gc - gm + pi*(sign(d1,(x - ga)) + (delta - d1)*sign(d1,(x - gb)) + sign(d1,x) + (d2 - deltas)*sign(d1,(x - gc)) - sign(d1,(x - gm))))))
end function 

! TODO: find true go -- gamma O
real(8) function v_g(x)
! ф. скорости V от гамма
use mod
real(8) x
!v_g = v_inf * dsin((x - go) * d5) / dabs(dsin((x - ga) * d5)) * (dabs(dsin((x - gc) * d5) / dsin((x - gb) * d5))) ** (delta - d1)
v_g = -v_inf * dsin((x - go) * d5) / dabs(dsin((x - ga) * d5)) * (dabs(dsin((x - gc) * d5) / dsin((x - gb) * d5))) ** (delta - d1)
end function

real(8) function s_g(a, b)
! дуговая абсцисса s(gamma)
! a, b -- соответственно нижний и верхний предел интегрирования
use mod
integer(4) i, nintv
real(8) a, b, res, ds_dgamma, errest
external ds_dgamma
if (a > b) then
    print *, "ERROR: a > b"
endif
call DQDAGS(ds_dgamma, a, b, d0, eps, res, errest)
s_g = res
end function

! TODO: implement, how to 
subroutine v_s
! процедура нахождения распр. скорости по пластинке и запись в файл
! распределение скорости по пластинке и запись в файл
use mod
integer(4) i, n
real(8) v_g, s_g, tmp, step, gj
external v_g, s_g
n = 301
!запись в файл 
open (port,FILE = 'vs.dat') 
write(port,*) ' VARIABLES = "S", "V" '
!до точки M, обратная сторона точки B
step = 0.01d0
tmp = d0
i = 1
write(port,"(F12.5,' ', F12.5)") d0, v_g(d2*pi)
do while(tmp <= l_2) 
    gj = d2 * pi - i * step
    tmp = s_g(gj, d2*pi)
    write(port,"(F12.5,' ', F12.5)") tmp, v_g(gj)    
    i = i + 1
enddo
!TODO: поделить на зоны
!с точки C, до конца
write(port,"(F12.5,' ', F12.5)") (l_2 + l_3), v_g(gc)
do i = 2, n
    gj = gc - (i - 1) * gc / (n - 1)
    !gj = d2 * pi - (i - 1) * d2 * pi / (n - 1) 
    write(port,"(F12.5,' ', F12.5)") (l_2 + l_3 + s_g(gj, gc)), v_g(gj)
end do
close(port)
port = port + 1
end subroutine

! TODO: implement
function z(x)
real(8) x
complex(8) z
z = c1
end function 

! TODO implement
subroutine find_shape_of_plast
! записывает в файл форму пластинки
use mod
integer(4) i 
complex(8) z, t

end subroutine
