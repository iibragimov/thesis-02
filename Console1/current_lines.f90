!============<Current lines calculation>==================
    
subroutine save_circle(port_zz)
!��������������� ��������� ��� ���������� ���������� � ����
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
!��������������� ��������� ��� ���������� ��������� � ����
use mod
integer(4) i, port_z
character(8) zone_name
complex(8) z_circle, t 
zone_name = 'plastina'
write(port_z, "('ZONE T=""', A5, i2, '"", I=', i4, ', F=POINT')") zone_name, 34, 3
!����� �������� ��� �������
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_a))), dimag(z(cdexp(ii * g_a)))
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_b))), dimag(z(cdexp(ii * g_b)))
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(c1)), dimag(z(c1))
t = l_2 * cdexp(ii * (d1 - delta) * pi)
write(port_z, "(F9.5, ' ', F9.5)") dreal(t), dimag(t)
write(port_z, "(F9.5, ' ', F9.5)") d0, d0
write(port_z, "(F9.5, ' ', F9.5)") -l_1, 0
end subroutine 

subroutine save_line(current_port, k, array, zone_name, zone_name_n)
!������ � ���� ����� ���� �� �����
!k -- ���-�� ����� ����� ����
!zl -- zl(0:nmax) ����������� ������ ����� ����
!zone_name -- �������� ����
!zone_name_n -- ����� �����
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
!���������� ������� dz_d\zeta, �������������� �������
!zz -- \zeta ����������� ����������
!�������� �� ������, ���� �� ����������
!��������� � ������ ������ ������
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
!���������� ������� dw_d\zeta
!zz -- \zeta ����������� ����������
!��������, ���� �� ����� ������� ���� ���������
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
!������� ��������� ��� ���������� ����� ����
!z -- ����� � ��������� z
!zz -- ����� � ��������� zeta
use mod
complex(8) z,zz
logical lines_test_stop
lines_test_stop = (cdabs(z) > 20.0d0) .OR. (cdabs(zz) > 20.0d0)
end function

function lines_test_stop_flow(zz, z)
use mod
!������� ��������� ��� ���������� ����� ����
!z -- ����� � ��������� z
!zz -- ����� � ��������� zeta
use mod
complex(8) z,zz
logical lines_test_stop
lines_test_stop = (dreal(z) > d0) .OR. (dreal(zz)) > d0
end function

function z_circle(g_up, g_down)
! z(\zeta) -- ����������� �� ������� �����������-������ �� ����������
! g_up -- �����, ������� ������
! g_down -- �����, ������ ������
use mod
integer(4) i
real(8) shag, g_up, g_down
complex(8) tn(nk), fn(nk) ! tn -- ����� ��������������, fn -- �������
complex(8) z_circle, dz_dzeta, cintz       
shag = (g_up - g_down) / (nk - 1)
do i = 1, nk
    tn(i) = cdexp(ii * (g_down + shag * (i - 1)))
    fn(i) = dz_dzeta(tn(i))
end do
z_circle = cintz(nk, tn, fn)
end function 

function z_plane(zz_up, zz_down)
! z(\zeta) -- ����������� �� ������� �����������-������ ��������� zeta �� z
! zz_up -- �������� ����� ��������������
! zz_down --  ��������� ����� ��������������
use mod
integer(4) i
complex(8) tn(nk), fn(nk) ! tn -- ����� ��������������, fn -- �������
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
! ���������� � ���� ����� ���������
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
!�������� ��������� (x, y, z)
use mod
external dz_dzeta, dw_dzeta, lines_test_stop, lines_test_stop_flow
logical lines_test_stop, lines_test_stop_flow
character(8) zone_name
integer(4) i, j, k, zone_name_n, number_of_lines, port_z, port_zz
real(8) dir0, dirz0, dl, shag, dk
complex(8) z0, zz0, zl(0:nmax), zlz(0:nmax), z_circle, z_plane, dw_dzeta, dz_dzeta, temp_zl_k
complex(8) zc, zm, tmp_z
!�������� ����������
!dw_dzeta -- dw_d\zeta �������
!dz_dzeta -- dz_d\zeta �������
!lines_test_stop -- ������� (�������) ��������� ������������� ��������, ������ ������� ����� ����
!number_of_lines -- ���������� ����� ����
!zz0 -- ��������� ����� ���������� ����� ���� zeta
!z0 -- ��������� ����� ���������� ����� ���� � ���������� z
!dir0 -- ��������� ����������� �������� ����� ���� � ���������� zeta dir0 = -zarg(dw_dzeta())
!dirz0 -- ��������� ����������� � ���������� z dir0 = -zarg(dz_dzeta())
!kdir - 1 - �� ������, -1 - ������ ������   
!upper_bound -- ������� ������� ����� ���� � ��������� zeta
!lower_bound -- ������ ������� ����� ���� � ��������� zeta
!zl - �������� ������ ����� � ��������� zeta
!zlz - �������� ������ ����� � ��������� z
!zone_name_n -- ������ ����� ����, (���� � ����� �� ������)
!port_z -- ������ ��� ������ � ���� ����� ���� � ��������� z
!port_zz -- ������ ��� ������ � ���� ����� ���� � ��������� zeta
!zc -- z � ����� � �� ���������
!zm -- z � ����� M �� ���������
!tmp_z -- ����� ��������� � ������� zlz � ������ �����
!���������� ����������/�������� ������
!call find_const !�������� �������� ������ � ������ ��� ���� ���������,����� �� ����� ����� ���������� xyz

!������������� ���������� ��� ���������� ����� ����
call init_lines_const

port_z = port
port_zz = port + 1
open(port_z, FILE='data/current_lines_z.dat')
open(port_zz, FILE='data/current_lines_zeta.dat')

!TODO:
!����� ����� ������ ����� ��, ������� ��������� � ��������� � �������������.
!����� �� ����� ����� ����� ����� ����������� ������ z0 � ����� � � ����������� ��� ��������������
!����� ������� ��������� � ��������������� ������� ��������������
!����� ������������ find_lines2

!1. �������� �� ����� ��������
zz0 = cdexp(ii * go)        !��������� zeta
z0 = -c1 + z_circle(go, ga) !��������� z
dir0 = go
dirz0 = 3/d2*pi
call find_line2(zz0, z0, dir0, dirz0, -1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
zone_name = 'first'
call save_line(port_z, k, zlz, zone_name, 0)
call save_line(port_zz, k, zl, zone_name, 0)

!2. �������� ����
temp_zl_k = zl(k)
number_of_lines = 20
shag = 0.02d0
tmp_z = zlz(k)

!��� ��������������� ������
do zone_name_n = 1, number_of_lines + 5
    zz0 = temp_zl_k * cdexp(-ii * shag * zone_name_n)   ! � ��������� zeta
    z0 = tmp_z + z_plane(zz0, temp_zl_k)                        ! � ��������� z
    dir0 = -zarg(dw_dzeta(zz0))
    dirz0 = -zarg(dz_dzeta(z0))
    call find_line2(zz0, z0, dir0, dirz0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
    zone_name = char(47 + zone_name_n)
    call save_line(port_z, k, zlz, zone_name, zone_name_n)
    call save_line(port_zz, k, zl, zone_name, zone_name_n)
end do
!��� ��������������� ������
do zone_name_n = 1, number_of_lines
    zz0 = temp_zl_k * cdexp(ii * shag * zone_name_n)    ! � ��������� zeta
    z0 = tmp_z + z_plane(zz0, temp_zl_k)                        ! � ��������� z
    dir0 = -zarg(dw_dzeta(zz0))
    dirz0 = -zarg(dz_dzeta(z0))
    call find_line2(zz0, z0, dir0, dirz0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
    zone_name = char(47 + zone_name_n)
    call save_line(port_z, k, zlz, zone_name, zone_name_n)
    call save_line(port_zz, k, zl, zone_name, zone_name_n)
end do

!3. �� ����� � ����� ������� 
zz0 = cdexp(ii * gc)                ! ��. zeta
zc = -c1 + z_circle(go, ga) + z_circle(gc+eps, go) 
z0 =  zc !��������� z
dir0 = gc
dirz0 = -pi * d5 * (delta - d1)    ! ��. z ����������� ���� � �. � �� ������ ����� (�����������)
call find_line2(zz0, z0, dir0, dirz0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
zone_name = 'razdel'
call save_line(port_z, k, zlz, zone_name, zone_name_n)
call save_line(port_zz, k, zl, zone_name, zone_name_n)

temp_zl_k = zl(k)
tmp_z = zlz(k)

!4. �� ������ ������
zz0 = c1
z0 = l_2 * cdexp(ii * (d1 - delta) * pi)
dir0 = d0 !-zarg(dw_dzeta(zz0))
dirz0 = -zarg(dz_dzeta(z0))
call find_line2(zz0, z0, dir0, dirz0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop)
zone_name = 'sxod'
call save_line(port_z, k, zlz, zone_name, zone_name_n)
call save_line(port_zz, k, zl, zone_name, zone_name_n)

!5. ����� ����� �� ������ ������ � �� ����� � ������ ������� ����
!��������� �������� �� ����� ������ ������ �������� ���� � ���������� � �������� ����������� ������ ������
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

!6. �������� � ����������
call save_plastin(port_z)
call save_circle(port_zz)

close(port_z)
close(port_zz)
port = port + 2
end subroutine
