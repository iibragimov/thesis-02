!============<Current lines calculation>==================
    
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
complex(8) z_circle, t 
zone_name = 'plastina'
write(port_z, "('ZONE T=""', A5, i2, '"", I=', i4, ', F=POINT')") zone_name, 34, 3
!нужно подумать как строить
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_a))), dimag(z(cdexp(ii * g_a)))
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_b))), dimag(z(cdexp(ii * g_b)))
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(c1)), dimag(z(c1))
t = l_2 * cdexp(ii * (d1 - delta) * pi)
write(port_z, "(F9.5, ' ', F9.5)") dreal(t), dimag(t)
write(port_z, "(F9.5, ' ', F9.5)") d0, d0
write(port_z, "(F9.5, ' ', F9.5)") -l_1, 0
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
dz_dzeta = (C * zz_a * zz_1 * zz_c * zz_c) * (zz_b / zz_c)**(delta) / (zz_b * zz_m * zz * zz)
!dz_dzeta = C * zz_a * zz_1 * zz_b**(delta - d1) * zz_c**(d2 - delta) / (zz_m * zz * zz)
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

function lines_test_stop_flow(zz, z)
use mod
!условие остановки при построении линии тока
!z -- точка в плоскости z
!zz -- точка в плоскости zeta
use mod
complex(8) z,zz
logical lines_test_stop
lines_test_stop = (dreal(z) > d0) .OR. (dreal(zz)) > d0
end function

function z_circle(g_up, g_down)
! z(\zeta) -- отображение по формуле Кристоффеля-Шварца на окружности
! g_up -- гамма, верхний предел
! g_down -- гамма, нижний предел
use mod
integer(4) i
real(8) shag, g_up, g_down
complex(8) tn(nk), fn(nk) ! tn -- линия интегрирования, fn -- функция
complex(8) z_circle, dz_dzeta, cintz       
shag = (g_up - g_down) / (nk - 1)
do i = 1, nk
    tn(i) = cdexp(ii * (g_down + shag * (i - 1)))
    fn(i) = dz_dzeta(tn(i))
end do
z_circle = cintz(nk, tn, fn)
end function 

function z_plane(zz_up, zz_down)
! z(\zeta) -- отображение по формуле Кристоффеля-Шварца плоскость zeta на z
! zz_up -- конечная точка интегрирования
! zz_down --  начальная точка интегрирования
use mod
integer(4) i
complex(8) tn(nk), fn(nk) ! tn -- линия интегрирования, fn -- функция
complex(8) z_plane, zz_up, zz_down, shag, dz_dzeta, cintz
shag = (zz_up - zz_down) / (nk - 1)
do i = 1, nk
    tn(i) = zz_down + shag * (i - 1)
    fn(i) = dz_dzeta(tn(i))
end do
z_plane = cintz(nk, tn, fn)
end function

! TODO implement or delete
subroutine find_shape_of_plast
! записывает в файл форму пластинки
use mod
integer(4) i 
real(8) g, up, down, shag
complex(8) z_circle, tmp
tmp = c0
shag = 0.050d0
down = gb
up = down + shag 
open(port, FILE='data/plastina.dat')

do i = 1, nk
    down = up
    up = up + shag
    tmp =  tmp + z_circle(up, down)
    write(port, "(F9.5, ' ', F9.5)") dreal(tmp), dimag(tmp)
end do

close(port)
print *, z_circle(ga, gb+eps)
port = port + 1
end subroutine

subroutine current_lines
!добавить параметры (x, y, z)
use mod
external dz_dzeta, dw_dzeta, lines_test_stop, lines_test_stop_flow
logical lines_test_stop, lines_test_stop_flow
character(8) zone_name
integer(4) i, j, k, zone_name_n, number_of_lines, port_z, port_zz
real(8) dir0, dirz0, dl, shag, dk
complex(8) z0, zz0, zl(0:nmax), zlz(0:nmax), z_circle, z_plane, dw_dzeta, dz_dzeta, temp_zl_k
complex(8) zc, zm, tmp_z
!Описание переменных
!dw_dzeta -- dw_d\zeta функция
!dz_dzeta -- dz_d\zeta функция
!lines_test_stop -- функции (условие) остановки итерационного процесса, задают границы линий тока
!number_of_lines -- количество линий тока
!zz0 -- начальная точка построения линии тока zeta
!z0 -- начальная точка построения линии тока в плосткости z
!dir0 -- начальное направление скорости линии тока в плосткости zeta dir0 = -zarg(dw_dzeta())
!dirz0 -- начальное направление в плосткости z dir0 = -zarg(dz_dzeta())
!kdir - 1 - по потоку, -1 - против потока   
!upper_bound -- верхняя граница линий тока в плоскости zeta
!lower_bound -- нижняя граница линий тока в плоскости zeta
!zl - выходной массив точек в плоскости zeta
!zlz - выходной массив точек в плоскости z
!zone_name_n -- индекс линий тока, (чтоб в файле не путать)
!port_z -- индекс для записи в файл линий тока в плоскости z
!port_zz -- индекс для записи в файл линий тока в плоскости zeta
!zc -- z в точке С на пластинке
!zm -- z в точке M на пластинке
!tmp_z -- точке последняя в массиве zlz в первой линии
!нахождение параметров/констант задачи
!call find_const !возможно придется убрать и искать вне этой процедуры,тогда не нужно будет передавать xyz

!инициализация параметров для нахождения линий тока
call init_lines_const

port_z = port
port_zz = port + 1
open(port_z, FILE='data/current_lines_z.dat')
open(port_zz, FILE='data/current_lines_zeta.dat')

!TODO:
!Нужно взять первую линию ту, которая врезается в пластинку и раздваивается.
!чтобы ее найти нужно найти точку раветвления потока z0 в точке О и направление для интегрирования
!Потом вызвать процедуру в противоположную сторону интегрирования
!нужно использовать find_lines2

!1. проинтег от точки разветвл
zz0 = cdexp(ii * go)        !плоскость zeta
z0 = -c1 + z_circle(go, ga) !плоскость z
dir0 = go
dirz0 = 3/d2*pi
call find_line2(zz0, z0, dir0, dirz0, -1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
zone_name = 'first'
call save_line(port_z, k, zlz, zone_name, 0)
call save_line(port_zz, k, zl, zone_name, 0)

!2. основной цикл
temp_zl_k = zl(k)
number_of_lines = 20
shag = 0.02d0
tmp_z = zlz(k)

!над разветлвяющейся линией
do zone_name_n = 1, number_of_lines + 5
    zz0 = temp_zl_k * cdexp(-ii * shag * zone_name_n)   ! в плоскости zeta
    z0 = tmp_z + z_plane(zz0, temp_zl_k)                        ! в плоскости z
    dir0 = -zarg(dw_dzeta(zz0))
    dirz0 = -zarg(dz_dzeta(z0))
    call find_line2(zz0, z0, dir0, dirz0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
    zone_name = char(47 + zone_name_n)
    call save_line(port_z, k, zlz, zone_name, zone_name_n)
    call save_line(port_zz, k, zl, zone_name, zone_name_n)
end do
!под разветлвяющейся линией
do zone_name_n = 1, number_of_lines
    zz0 = temp_zl_k * cdexp(ii * shag * zone_name_n)    ! в плоскости zeta
    z0 = tmp_z + z_plane(zz0, temp_zl_k)                        ! в плоскости z
    dir0 = -zarg(dw_dzeta(zz0))
    dirz0 = -zarg(dz_dzeta(z0))
    call find_line2(zz0, z0, dir0, dirz0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
    zone_name = char(47 + zone_name_n)
    call save_line(port_z, k, zlz, zone_name, zone_name_n)
    call save_line(port_zz, k, zl, zone_name, zone_name_n)
end do

!3. от точки С линия раздела 
zz0 = cdexp(ii * gc)                ! пл. zeta
zc = -c1 + z_circle(go, ga) + z_circle(gc+eps, go) 
z0 =  zc !плоскость z
dir0 = gc
dirz0 = -pi * d5 * (delta - d1)    ! пл. z биссектриса угла в т. С со знаком минус (направление)
call find_line2(zz0, z0, dir0, dirz0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
zone_name = 'razdel'
call save_line(port_z, k, zlz, zone_name, zone_name_n)
call save_line(port_zz, k, zl, zone_name, zone_name_n)

temp_zl_k = zl(k)
tmp_z = zlz(k)

!4. от задней кромки
zz0 = c1
z0 = l_2 * cdexp(ii * (d1 - delta) * pi)
dir0 = d0 !-zarg(dw_dzeta(zz0))
dirz0 = -zarg(dz_dzeta(z0))
call find_line2(zz0, z0, dir0, dirz0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
zone_name = 'sxod'
call save_line(port_z, k, zlz, zone_name, zone_name_n)
call save_line(port_zz, k, zl, zone_name, zone_name_n)

!5. выдув струи от задней кромки и от точки С линией раздела сред
!вычисляем значения от линии задней кромки отступаю шаги и интегрируя в обратном направлении против потока
number_of_lines = 3
shag = (dimag(zl(k)) - dimag(temp_zl_k)) / number_of_lines
temp_zl_k = zl(k)
tmp_z = zlz(k)
do i = 1, number_of_lines + 2
    zz0 = temp_zl_k - ii*shag*i
    z0 = tmp_z + z_plane(zz0, temp_zl_k)
    dir0 = -zarg(dw_dzeta(zz0))
    dirz0 = -zarg(dz_dzeta(z0))
    call find_line2(zz0, z0, dir0, dirz0, -1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop_flow)
    do j = 0, k
        if (dimag(zlz(j)) > 0) then
            zlz(j) = dreal(zlz(j)) * c1
        end if
    end do
    zone_name = 'flow'
    call save_line(port_z, k, zlz, zone_name, zone_name_n)
    call save_line(port_zz, k, zl, zone_name, zone_name_n)
end do

!6. пластина и окружность
call save_plastin(port_z)
call save_circle(port_zz)

close(port_z)
close(port_zz)
port = port + 2
end subroutine
