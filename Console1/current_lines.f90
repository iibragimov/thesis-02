!============<Current lines calculation>==================

! TODO implement or delete
subroutine find_shape_of_plast
! записывает в файл форму пластинки
use mod
integer(4) i 
real(8) shag, g
complex(8) z, tmp
shag = (ga - gb) / (nmax - 1)
open(port, FILE='data/plastina.dat')
do i = 1, nmax
    g = gb + shag * (i - 1)
    tmp = z(cdexp(ii * g))
    write(port, "(F9.5, ' ', F9.5)") dreal(tmp), dimag(tmp)
end do
close(port)
port = port + 1
end subroutine
    
subroutine save_circle(port_zz)
!вспомогательная процедура для сохранения окружности в файл
use mod
integer(4) i, port_zz, k
character(8) zone_name 
real(8) gj
zone_name = 'circle'
k = 201
write(port_zz, "('ZONE T=""', A5, i2, '"", I=', i4, ', F=POINT')") zone_name, 35, 201
do i = 1, k
    gj = d2 * pi - (i - 1) * d2 * pi / (k - 1)
    write(port_zz, "(F9.5, ' ', F9.5)") dcos(gj), dsin(gj)
end do
end subroutine 

subroutine save_plastin(port_z)
!вспомогательная процедура для сохранения пластинки в файл
use mod
integer(4) i, port_z
character(8) zone_name
complex(8) z, t 
zone_name = 'plastina'
write(port_z, "('ZONE T=""', A5, i2, '"", I=', i4, ', F=POINT')") zone_name, 34, 3
!нужно подумать как строить
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_a))), dimag(z(cdexp(ii * g_a)))
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_b))), dimag(z(cdexp(ii * g_b)))
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(c1)), dimag(z(c1))
end subroutine 

subroutine save_line(current_port, k, array, zone_name, zone_name_n)
!запись в файл линии тока по зонам
!k -- кол-во точек линии тока
!zl -- zl(0:nmax) комплексный массив линии тока
!zone_name -- название зоны
!zone_name_n -- номер линии
use mod 
integer(4) i, k, zone_name_n, current_port
character(8) zone_name
complex(8) array(0:nmax)
write(current_port, "('ZONE T=""', A5, i2, '"", I=', i4, ', F=POINT')") zone_name, zone_name_n, k+1
do i = 0, k
    write(current_port, "(F9.5, ' ', F9.5)") dreal(array(i)), dimag(array(i))
enddo
end subroutine 

function dz_dzeta(zz)
!нахождение функции dz_d\zeta, поынтегральная функция
!zz -- \zeta комплексная переменная
!разделил на скобки, чтоб не запупаться
!разделить с учетом выбора ветвей
use mod
complex(8) zz, dz_dzeta
complex(8) C, zz_a, zz_b, zz_1, zz_m, zz_c
C = mod_C * cdexp(ii * arg_C)
zz_a = (zz - cdexp(ii * ga))
zz_b = (zz - cdexp(ii * gb))
zz_c = (zz - cdexp(ii * gc))
zz_m = (zz - cdexp(ii * gm))
zz_1 = (zz - c1)
!dz_dzeta = (C * zz_a * zz_1 * zz_c * zz_c) * (zz_b / zz_c)**(delta) / (zz_b * zz_m * zz * zz)
dz_dzeta = C * zz_a * zz_1 * zz_b**(delta - d1) * zz_c**(d2 - delta) / (zz_m * zz * zz)
end function

function dw_dzeta(zz)
!нахождение функции dw_d\zeta
!zz -- \zeta комплексная переменная
!разделил, чтоб не таким длинным было выражение
use mod
complex(8) zz, dw_dzeta
complex(8) u, zz_o, zz_c, zz_1, zz_m
u = u_inf * cdexp(-ii * beta)
zz_o = (zz - cdexp(ii * go))
zz_c = (zz - cdexp(ii * gc))
zz_m = (zz - cdexp(ii * gm))
zz_1 = (zz - c1) 
dw_dzeta = (u * zz_o * zz_c * zz_1) / (zz_m * zz * zz)
end function

function lines_test_stop(zz, z)
use mod
!условие остановки при построении линии тока
!z -- точка в плоскости z
!zz -- точка в плоскости zeta
use mod
complex(8) z,zz
logical lines_test_stop
lines_test_stop = (cdabs(z) > 20.0d0) .OR. (cdabs(zz) > 20.0d0)
end function

!TODO: implement
function z(zz)
! z(\zeta) -- отображение по формуле Кристоффеля-Шварца
! Берем интеграл от gb до то точки zz (zeta)
use mod
integer(4) i
real(8) shag
complex(8) t(nmax), f(nmax) ! t -- линия интегрирования, f -- функция
complex(8) z, zz, dz_dzeta, cintz
shag = (zarg(zz) - gb) / (nmax - 1)
do i = 1, nmax
    t(i) = cdexp(ii * (gb + shag * (i - 1)))
    f(i) = dz_dzeta(t(i))
end do
z = cintz(nmax, t, f)
end function 

subroutine current_lines
!добавить параметры (x, y, z)
use mod
external dz_dzeta, dw_dzeta, lines_test_stop
logical lines_test_stop
character(8) zone_name
integer(4) i, k, zone_name_n, number_of_lines, port_z, port_zz
real(8) dir0, dl, shag, dk
complex(8) z0, zz0, zl(0:nmax), zlz(0:nmax), z, dw_dzeta, dz_dzeta, temp_zl_k
!Описание переменных
!dw_dzeta -- dw_d\zeta функция
!dz_dzeta -- dz_d\zeta функция
!lines_test_stop -- функции (условие) остановки итерационного процесса, задают границы линий тока
!number_of_lines -- количество линий тока
!zz0 -- начальная точка построения линии тока zeta
!z0 -- начальная точка построения линии тока в плосткости z
!dir0 -- начальное направление скорости линии тока в плосткости zeta dir0 = -zarg(dw_dzeta())
!kdir - 1 - по потоку, -1 - против потока   
!upper_bound -- верхняя граница линий тока в плоскости zeta
!lower_bound -- нижняя граница линий тока в плоскости zeta
!zl - выходной массив точек в плоскости zeta
!zlz - выходной массив точек в плоскости z
!zone_name_n -- индекс линий тока, (чтоб в файле не путать)
!port_z -- индекс для записи в файл линий тока в плоскости z
!port_zz -- индекс для записи в файл линий тока в плоскости zeta

!нахождение параметров/констант задачи
!call find_const !возможно придется убрать и искать вне этой процедуры,тогда не нужно будет передавать xyz

!инициализация параметров для нахождения линий тока
call init_lines_const

!TODO:
!Нужно взять первую линию ту, которая врезается в пластинку и раздваивается.
!чтобы ее найти нужно найти точку раветвления потока z0 в точке О и направление для интегрирования
!Потом вызвать процедуру в противоположную сторону интегрирования
!нужно использовать find_lines2

port_z = port
port_zz = port + 1
open(port_z, FILE='data/current_lines_z.dat')
open(port_zz, FILE='data/current_lines_zeta.dat')

!1. проинтег от точки разветвл
!zz0 =
!2. основной цикл
!3. от задней кромки

call save_plastin(port_z)
call save_circle(port_zz)

close(port_z)
close(port_zz)
port = port + 2
end subroutine
