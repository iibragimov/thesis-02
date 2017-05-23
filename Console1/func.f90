!function sred(a,b,fa,fb,c)
!! среднее
!real(8) sred,a,b,fa,fb,c
!sred=fa+(fb-fa)/(b-a)*(c-a)
!end

function csred(a,b,fa,fb,c)
! среднее
real(8) a,b,c
complex(8) csred,fa,fb
csred=fa+(fb-fa)/(b-a)*(c-a)
end

function scal_p(a,b)
!скал€рное произведение
real(8) scal_p
complex(8) a,b
scal_p=dreal(dconjg(a)*b)
end

function vect_p(a,b)
!векторное произведение
real(8) vect_p
complex(8) a,b
vect_p=dimag(dconjg(a)*b)
end

function get_s_section(z,z1,z2)
!вычисление дуговой абсциссы проекции точки z на отрезок (z1,z2)
!s=0 - проекци€ совпадает с z1
!s=1 - проекци€ совпадает с z2
complex(8) z,z1,z2,dz21,dz1
real(8) get_s_section
dz21=z2-z1
dz1=z-z1
get_s_section=dreal(dconjg(dz21)*dz1)/cdabs(dz21)**2
end

subroutine sort_mass(n,a)
!сотрировка массива
integer(4) i,j,n,jmin
real(8) a(n),tmin
do i=1,n-1
  tmin=a(i)
  jmin=i
  do j=i+1,n
    if (a(j)<tmin) then
	  tmin=a(j)
	  jmin=j
	endif
  enddo
  a(jmin)=a(i)
  a(i)=tmin
enddo
end

subroutine angles_razryv(nj1,alpha)
!удаление разрыва при скачках функций на 2\pi
use gen_mod
integer(4) i,nj1
real(8) dd,alpha(nj1),dpi
  dpi=pi !*1.5d0
  dd=d0
  do i=2,nj1
    alpha(i)=alpha(i)+dd
    if ((alpha(i)-alpha(i-1))>dpi) then 
      dd=dd-pi2
	  alpha(i)=alpha(i)-pi2
    endif
    if ((alpha(i)-alpha(i-1))<-dpi) then 
      dd=dd+pi2
	  alpha(i)=alpha(i)+pi2
    endif
  enddo
end

subroutine spline_razryv(xx,yy,bb,cc,nr,nmax,razryv)
!построение сплайна от функции с разрывной производной
use gen_mod
integer(4) i,j,k,nr,nmax,razryv(nmax)
real(8) xx(nmax),yy(nmax),bb(nmax),cc(4,nmax)
do i=1,nr-1
  j=razryv(i)
  k=razryv(i+1)
  call dcsakm(k-j+1,xx(j:k),yy(j:k),bb(j:k),cc(:,j:k))
enddo
end

subroutine spline_linear(nj1,xx,yy,bb,cc)
!линейный сплайн
use gen_mod
integer(4) i,nj1,j,jmin
real(8) xx(nj1),yy(nj1),bb(nj1),cc(4,nj1),yy1(nj1)
real(8) tmin,fmin
bb=xx
yy1=yy
do i=1,nj1-1
  tmin=bb(i)
  fmin=yy1(i)
  jmin=i
  do j=i+1,nj1
    if (bb(j)<tmin) then
	  tmin=bb(j)
	  fmin=yy1(j)
	  jmin=j
	endif
  enddo
  bb(jmin)=bb(i)
  bb(i)=tmin
  yy1(jmin)=yy1(i)
  yy1(i)=fmin
enddo
  
  cc=d0
  do i=1,nj1-1
    cc(1,i)=yy1(i) 
    cc(2,i)=(yy1(i+1)-yy1(i))/(bb(i+1)-bb(i))
  enddo
  cc(1,nj1)=yy1(nj1)
  cc(2,nj1)=cc(2,nj1-1)
end

subroutine findline(nt,tj,tt,x0,y0,lx,ly)
!определение формы линии по углу наклона касательной
use gen_mod
integer(4) nt       !число точек
real(8) tj(nt)      !дугова€ абсцисса 
real(8) tt(nt)      !угол наклона касательной к кривой
real(8) x0,y0       !начальна€ точка
real(8) lx(nt),ly(nt) !координаты точек кривой
complex(8) f(nt)

f=c1
call findlinez(nt,tj,tt,f,x0,y0,lx,ly)
end

subroutine findlinez(nt,tj,tt,f,x0,y0,lx,ly)
!определение формы линии интегрированием уравнени€ dz = f(\zeta)d\zeta
use gen_mod
integer(4) nt       !число точек
real(8) tj(nt)      !дугова€ абсцисса 
real(8) tt(nt)      !угол наклона касательной к кривой
complex(8) f(nt)    !функци€ f(\zeta)
real(8) x0,y0       !начальна€ точка
real(8) lx(nt),ly(nt) !координаты точек кривой
integer(4) i
complex(8) fc(nt)

do i=1,nt
  fc(i)=f(i)*cdexp(ii*tt(i))
enddo
call findlinez1(nt,tj,fc,x0,y0,lx,ly)
end

subroutine findlinez1(nt,tj,f,x0,y0,lx,ly)
!определение формы линии интегрированием уравнени€ dz = f(\zeta)d\sigma
!d\sigma=|d\zeta|
use gen_mod
integer(4) nt       !число точек
real(8) tj(nt)      !дугова€ абсцисса 
complex(8) f(nt)    !функци€ f(\zeta)
real(8) x0,y0       !начальна€ точка
real(8) lx(nt),ly(nt) !координаты точек кривой
integer(4) i
real(8) ff(nt)
real(8) bb(nt),cc(4,nt)

ff=dreal(f)
call dcsakm(nt,tj,ff,bb,cc)
lx(1)=x0
do i=2,nt
  lx(i)=lx(i-1)+dcsitg(tj(i-1),tj(i),nt-1,bb,cc)
enddo
ff=dimag(f)
call dcsakm(nt,tj,ff,bb,cc)
ly(1)=y0
do i=2,nt
  ly(i)=ly(i-1)+dcsitg(tj(i-1),tj(i),nt-1,bb,cc)
enddo
end

function cint(nt,tj,f,t0,t1)
!вычисление комплексного интеграла \int_t0^t1 f(t) dt
use gen_mod
complex(8) cint
integer(4) nt       !число точек
real(8) tj(nt)      !дугова€ абсцисса 
complex(8) f(nt)    !функци€ f(t)
real(8) t0,t1       !пределы интегрировани€ 
real(8) d_int
cint=c1*d_int(nt,tj,dreal(f),t0,t1)+ii*d_int(nt,tj,dimag(f),t0,t1)
end

function d_int(nt,tj,f,t0,t1)
!вычисление вещественного интеграла \int_t0^t1 f(t) dt
use gen_mod
real(8) d_int
integer(4) nt       !число точек
real(8) tj(nt)      !дугова€ абсцисса 
real(8) f(nt)    !функци€ f(t)
real(8) t0,t1       !пределы интегрировани€ 
real(8) bb(nt),cc(4,nt)
call dcsakm(nt,tj,f,bb,cc)
d_int=dcsitg(t0,t1,nt-1,bb,cc)
end

function cintz(nt,t,f)
!вычисление комплексного интеграла \int_t0^t1 f(t) dt
use gen_mod
complex(8) cintz
integer(4) nt       !число точек
complex(8) t(nt)    !координаты линии интегрировани€
complex(8) f(nt)    !функци€ f(t)
complex(8) cint
real(8) tj(nt)      !дугова€ абсцисса 
complex(8) f1(nt)
real(8) xx(nt),yy(nt),bb(nt),cc(4,nt)
integer(4) i
tj(1)=d0
do i=2,nt
  tj(i)=tj(i-1)+cdabs(t(i)-t(i-1))
enddo
do i=1,nt
  xx(i)=dreal(t(i))
  yy(i)=dimag(t(i))
enddo
call dcsakm(nt,tj,xx,bb,cc)
do i=1,nt
  xx(i)=dcsder(1,tj(i),nt-1,bb,cc)
enddo
call dcsakm(nt,tj,yy,bb,cc)
do i=1,nt
  yy(i)=dcsder(1,tj(i),nt-1,bb,cc)
enddo
do i=1,nt
  f1(i)=f(i)*(c1*xx(i)+ii*yy(i))
enddo
cintz=cint(nt,tj,f1,tj(1),tj(nt))
end

