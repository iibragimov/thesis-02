module gen_mod
! подключение IMSL
include 'link_fnl_shared.h'
use LSLRG_INT
! константы
real(8),PARAMETER :: d0=0.0d0
real(8),PARAMETER :: d1=1.0d0
real(8),PARAMETER :: d2=2.0d0
real(8),PARAMETER :: d5=0.50d0
real(8),PARAMETER :: pi=3.14159265358979323d0
real(8),PARAMETER :: pi2=pi*d2
real(8),PARAMETER :: pi5=pi*d5
complex(8),PARAMETER :: ii=(0.0d0,1.0d0)
complex(8),PARAMETER :: c0=(0.0d0,0.0d0)
complex(8),PARAMETER :: c1=(1.0d0,0.0d0)
complex(8),PARAMETER :: c2=(2.0d0,0.0d0)
complex(8),PARAMETER :: c3=(3.0d0,0.0d0)
real(8),PARAMETER :: eps = 0.000010d0          ! epsulon besk. malaya
end module gen_mod
