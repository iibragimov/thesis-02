!============<Additional functions, not important (ненужные функции)>==================

!вроде не нужна
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

!вроде тоже не нужна
function zzeta(x, y)
! ф-ция написана для сокращения кода при вычислении dz(g)
! разложение (zeta - zeta0) - с учетом выбора ветвей
use mod
real(8) x, y
complex(8) zzeta
zzeta = cdexp(ii * (x + y + pi * dsign(d1, x - y)) / d2) * dabs(d2 * dsin((x - y) / d2))
!zzeta = cdexp(ii * (x + y + pi) / d2) * dabs(d2 * dsin((x - y) / d2))
end function 

!вроде тоже не нужна
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

!вроде тоже не нужна
function sgu(x)
! сгущает точки к середине
use mod
real(8) x, sgu
sgu = (x - d5) * (x - d5) * (x - d5) / (d5 * d5) + d5
!sgu=d1 - ( atan( tan(( (1-x)*(1-x) )*pi/d2) ) )*d2/pi
end function    

!вроде тоже не нужна
function dz_dzeta_okr(x)
! поынтегральная ф. от вещественной переменной, скорее всего для вычисления на окружности
! упрощал эту функцию, поэтому возможны ошибки 09.04.2017
use mod
real(4) x
complex(8) C, dz_dzeta_okr
dz_dzeta_okr = mod_C * dabs(d2*dsin(d5*(x-ga))) * (dabs(d2*dsin(d5*(x-gb))) ** (delta - d1)) * dabs(d2*dsin(d5*x)) * (dabs(d2*dsin(d5*(x-gc))) ** (d2 - delta)) / dabs(d2*dsin(d5*(x-gm))) * cdexp(ii*(arg_C - g + d5*(ga + (delta - d1)*gb + (d2 - delta)*gc - gm + pi*(sign(d1,(x - ga)) + (delta - d1)*sign(d1,(x - gb)) + sign(d1,x) + (d2 - deltas)*sign(d1,(x - gc)) - sign(d1,(x - gm))))))
end function 