!============<Params of the task (velocity, density, coef of lifting force)>==================

! TODO: find true go -- gamma O
real(8) function v_g(x)
! ф. скорости V от гамма
use mod
real(8) x
!v_g = v_inf * dsin((x - go) * d5) / dabs(dsin((x - ga) * d5)) * (dabs(dsin((x - gc) * d5) / dsin((x - gb) * d5))) ** (delta - d1)
v_g = -v_inf * dsin((x - go) * d5) / dabs(dsin((x - ga) * d5)) * (dabs(dsin((x - gc) * d5) / dsin((x - gb) * d5))) ** (delta - d1)
end function

!TODO: оптимизировать вычисление интеграла (складывая по чуть-чуть)
real(8) function s_g(a, b)
! дуговая абсцисса s(gamma)
! a, b -- соответственно нижний и верхний предел интегрирования
use mod
integer(4) i, nintv
real(8) a, b, res, ds_dgamma, errest
external ds_dgamma
if (a > b) then
    print *, "ERROR in function s_g(a, b): a > b"
endif
call DQDAGS(ds_dgamma, a, b, d0, eps, res, errest)
s_g = res
end function

! TODO: учесть знаки, чтобы правильно вычислялось 
subroutine v_s
! процедура нахождения распр. скорости по пластинке и запись в файл
! распределение скорости по пластинке и запись в файл
use mod
integer(4) i, n
real(8) v_g, s_g, tmp, step, gj
external v_g, s_g
n = 301
!запись в файл 
open (port,FILE = 'data/vs.dat') 
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
! TODO: поделить на зоны
! с точки C, до конца
write(port,"(F12.5,' ', F12.5)") (l_2 + l_3), v_g(gc)
do i = 2, n
    gj = gc - (i - 1) * gc / (n - 1)
    !gj = d2 * pi - (i - 1) * d2 * pi / (n - 1) 
    write(port,"(F12.5,' ', F12.5)") (l_2 + l_3 + s_g(gj, gc)), v_g(gj)
end do
close(port)
port = port + 1
end subroutine
