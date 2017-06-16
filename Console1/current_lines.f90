!============<Current lines calculation>==================

! TODO implement or delete
subroutine find_shape_of_plast
! ���������� � ���� ����� ���������
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
complex(8) z, t 
zone_name = 'plastina'
write(port_z, "('ZONE T=""', A5, i2, '"", I=', i4, ', F=POINT')") zone_name, 34, 3
!����� �������� ��� �������
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_a))), dimag(z(cdexp(ii * g_a)))
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_b))), dimag(z(cdexp(ii * g_b)))
!write(port_z, "(F9.5, ' ', F9.5)") dreal(z(c1)), dimag(z(c1))
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
!dz_dzeta = (C * zz_a * zz_1 * zz_c * zz_c) * (zz_b / zz_c)**(delta) / (zz_b * zz_m * zz * zz)
dz_dzeta = C * zz_a * zz_1 * zz_b**(delta - d1) * zz_c**(d2 - delta) / (zz_m * zz * zz)
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

!TODO: implement
function z(zz)
! z(\zeta) -- ����������� �� ������� �����������-������
! ����� �������� �� gb �� �� ����� zz (zeta)
use mod
integer(4) i
real(8) shag
complex(8) t(nmax), f(nmax) ! t -- ����� ��������������, f -- �������
complex(8) z, zz, dz_dzeta, cintz
shag = (zarg(zz) - gb) / (nmax - 1)
do i = 1, nmax
    t(i) = cdexp(ii * (gb + shag * (i - 1)))
    f(i) = dz_dzeta(t(i))
end do
z = cintz(nmax, t, f)
end function 

subroutine current_lines
!�������� ��������� (x, y, z)
use mod
external dz_dzeta, dw_dzeta, lines_test_stop
logical lines_test_stop
character(8) zone_name
integer(4) i, k, zone_name_n, number_of_lines, port_z, port_zz
real(8) dir0, dl, shag, dk
complex(8) z0, zz0, zl(0:nmax), zlz(0:nmax), z, dw_dzeta, dz_dzeta, temp_zl_k
!�������� ����������
!dw_dzeta -- dw_d\zeta �������
!dz_dzeta -- dz_d\zeta �������
!lines_test_stop -- ������� (�������) ��������� ������������� ��������, ������ ������� ����� ����
!number_of_lines -- ���������� ����� ����
!zz0 -- ��������� ����� ���������� ����� ���� zeta
!z0 -- ��������� ����� ���������� ����� ���� � ���������� z
!dir0 -- ��������� ����������� �������� ����� ���� � ���������� zeta dir0 = -zarg(dw_dzeta())
!kdir - 1 - �� ������, -1 - ������ ������   
!upper_bound -- ������� ������� ����� ���� � ��������� zeta
!lower_bound -- ������ ������� ����� ���� � ��������� zeta
!zl - �������� ������ ����� � ��������� zeta
!zlz - �������� ������ ����� � ��������� z
!zone_name_n -- ������ ����� ����, (���� � ����� �� ������)
!port_z -- ������ ��� ������ � ���� ����� ���� � ��������� z
!port_zz -- ������ ��� ������ � ���� ����� ���� � ��������� zeta

!���������� ����������/�������� ������
!call find_const !�������� �������� ������ � ������ ��� ���� ���������,����� �� ����� ����� ���������� xyz

!������������� ���������� ��� ���������� ����� ����
call init_lines_const

!TODO:
!����� ����� ������ ����� ��, ������� ��������� � ��������� � �������������.
!����� �� ����� ����� ����� ����� ����������� ������ z0 � ����� � � ����������� ��� ��������������
!����� ������� ��������� � ��������������� ������� ��������������
!����� ������������ find_lines2

port_z = port
port_zz = port + 1
open(port_z, FILE='data/current_lines_z.dat')
open(port_zz, FILE='data/current_lines_zeta.dat')

!1. �������� �� ����� ��������
!zz0 =
!2. �������� ����
!3. �� ������ ������

call save_plastin(port_z)
call save_circle(port_zz)

close(port_z)
close(port_zz)
port = port + 2
end subroutine
