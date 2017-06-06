module mod

! include IMSL and constants
use gen_mod

! port number for writing to file. 
! need to increment ++ after recording to file
integer(4) port 

! main parameters of problem
real(8),PARAMETER :: l_1 = 1.0d0                ! lenght of plate / длина пластинки
real(8),PARAMETER :: l_2 = 0.10d0               ! lenght of flap / длина закрылка
real(8),PARAMETER :: v_inf = 1.0d0              ! velocity in infinity / скорость на бесконечности
real(8),PARAMETER :: delta = 1.30d0             ! coef of the deflection angle of the flap / коэф угола отклонения закрылка delta = [1, 1.5], (delta * pi -- угол)
real(8),PARAMETER :: alpha = pi/6               ! angle of attack / угол атаки
real(8),PARAMETER :: h = 0.10d0                 ! width of channel / ширина канала
real(8),PARAMETER :: l_3=h/dsin(delta*pi-pi)    ! leght of gap / длина щели

integer(4),PARAMETER :: nmax = 201

! unknown parameters
real(8) beta, u_inf                             ! angle of attack and velocity in inf in parametric plane / угол атаки и скорость на бесконечности в параметрической плоскости
real(8) mod_C, arg_C                            ! module and argument of complex constsnt / модуль и аргумент комплексной константы
real(8) ga, gb, gc, gm, go                      ! unknown angles / неизвестные углы
real(8) Q, Gamma, Cy                            ! flow and circulation / расход и циркуляция
complex(8) tt(nmax)                             ! form of plate / форма пластинки

end module mod