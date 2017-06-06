module mod

! include IMSL and constants
use gen_mod

! port number for writing to file. 
! need to increment ++ after recording to file
integer(4) port 

! main parameters of problem
real(8),PARAMETER :: l_1 = 1.0d0                ! lenght of plate / ����� ���������
real(8),PARAMETER :: l_2 = 0.10d0               ! lenght of flap / ����� ��������
real(8),PARAMETER :: v_inf = 1.0d0              ! velocity in infinity / �������� �� �������������
real(8),PARAMETER :: delta = 1.30d0             ! coef of the deflection angle of the flap / ���� ����� ���������� �������� delta = [1, 1.5], (delta * pi -- ����)
real(8),PARAMETER :: alpha = pi/6               ! angle of attack / ���� �����
real(8),PARAMETER :: h = 0.10d0                 ! width of channel / ������ ������
real(8),PARAMETER :: l_3=h/dsin(delta*pi-pi)    ! leght of gap / ����� ����

integer(4),PARAMETER :: nmax = 201

! unknown parameters
real(8) beta, u_inf                             ! angle of attack and velocity in inf in parametric plane / ���� ����� � �������� �� ������������� � ��������������� ���������
real(8) mod_C, arg_C                            ! module and argument of complex constsnt / ������ � �������� ����������� ���������
real(8) ga, gb, gc, gm, go                      ! unknown angles / ����������� ����
real(8) Q, Gamma, Cy                            ! flow and circulation / ������ � ����������
complex(8) tt(nmax)                             ! form of plate / ����� ���������

end module mod