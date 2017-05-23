module mod
! ����������� �������� � IMSL
use gen_mod
! ����� ����� ��� ������ � ����. 
! �� ����� ��� ����� ������������� ��� ������ � ����, ��� ����� ����� ����� ����������
integer(4) port 
! �������� ��������� ������ 
real(8),PARAMETER :: l_1 = 1.0d0                ! ����� ���������
real(8),PARAMETER :: l_2 = 0.10d0               ! ����� ��������
real(8),PARAMETER :: v_inf = 1.0d0              ! �������� �� �������������
real(8),PARAMETER :: delta = 1.30d0             ! ���� ����� ���������� �������� delta = [1, 1.5], (delta * pi -- ����)
real(8),PARAMETER :: alpha = pi / 6             ! ���� �����
real(8),PARAMETER :: h = 0.10d0                 ! ������ ������
real(8),PARAMETER :: l_3=h/dsin(delta*pi-pi)         ! ����� ���� 

integer(4),PARAMETER :: nmax = 201

! ����������� ���������
real(8) beta, u_inf                             ! ���� ����� � �������� �� ������������� � ��������������� ���������
real(8) mod_C, arg_C                            ! ������ � �������� ����������� ��������� 
real(8) ga, gb, gc, gm, go                      ! ����������� ���� 
real(8) Q, Gamma, Cy                            ! ������ � ����������
complex(8) tt(nmax)                             ! ����� ���������

end module mod