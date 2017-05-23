module mod
! подключение констант и IMSL
use gen_mod
! номер порта для записи в файл. 
! он будет все время увеличиваться при записи в файл, тем самым порты будут различными
integer(4) port 
! основные параметры задачи 
real(8),PARAMETER :: l_1 = 1.0d0                ! длина пластинки
real(8),PARAMETER :: l_2 = 0.10d0               ! длина закрылка
real(8),PARAMETER :: v_inf = 1.0d0              ! скорость на бесконечности
real(8),PARAMETER :: delta = 1.30d0             ! коэф угола отклонения закрылка delta = [1, 1.5], (delta * pi -- угол)
real(8),PARAMETER :: alpha = pi / 6             ! угол атаки
real(8),PARAMETER :: h = 0.10d0                 ! ширина канала
real(8),PARAMETER :: l_3=h/dsin(delta*pi-pi)         ! длина щели 

integer(4),PARAMETER :: nmax = 201

! неизвестные параметры
real(8) beta, u_inf                             ! угол атаки и скорость на бесконечности в параметрической плоскости
real(8) mod_C, arg_C                            ! модуль и аргумент комплексной константы 
real(8) ga, gb, gc, gm, go                      ! неизвестный углы 
real(8) Q, Gamma, Cy                            ! Расход и Циркуляция
complex(8) tt(nmax)                             ! форма пластинки

end module mod