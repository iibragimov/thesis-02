!��������� ���������� ����� ����
!subroutine init_lines_const    !������������� ��������   
!subroutine find_d_direction(f_dwdz,z0,a0,b0,c)  !���������� ����������� ��������� �� ����������� ����� ����� ���� (������� ����������� �����, ���� ������� �������)
!subroutine find_line(z0,dir0,kdir,zl,k,nmax,f_dwdzeta)  !���������� ����� ���� � ��������� zeta
!subroutine find_line2(z0,zz0,dir0,dirz0,kdir,zl,zlz,k,nmax,f_dwdzeta,f_dzdzeta,f_test_stop)  !���������� ����� ���� � ��������� zeta � z

!!!����� ������������� ���������� line_mod.f90 � func.f90
    
subroutine init_lines_const
!������������� ��������
use lines_mod
lines_eps_d_dir=1.0d-6 
lines_kk=1.5d0
lines_eps_dir=0.05d0
lines_dsmax=0.1d0
lines_dsmin=0.0001d0
lines_rmax=5.0d0
lines_d_dir_dr=1.0d-5
end
    
subroutine find_d_direction(f_dwdz,z0,a0,b0,c)
!���������� ����������� ��������� �� ����������� ����� ����� ���� (������� ����������� �����, ���� ������� �������)
!f_dwdz ������� ������� ���������� ����������� �������� [complex(8) function f_dwdz(complex(8))]
!z0 - ���������� ����������� �����
!a0,b0 - �������� ����� � ������� ������� ����� (���� a=0 � b=0, �������� ���������� a=0 � b=pi)
!c - �������� �������� ���������� ���� 
use gen_mod
use lines_mod
external f_dwdz
complex(8) f_dwdz,z0
real(8) a,b,fa,fb,c,fc,f_d_dir,a0,b0

a=a0
b=b0
if (a==d0.and.b==d0) b=pi

fa=f_d_dir(a,f_dwdz,z0)
fb=f_d_dir(b,f_dwdz,z0)
do while (b-a>lines_eps_d_dir)
    c=(a+b)*d5
    fc=f_d_dir(c,f_dwdz,z0)
    if (fa*fc>0) then
        a=c
        fa=fc
    else
        b=c
        fb=fc
    endif
enddo    
end

function f_d_dir(g,f_dwdz,z0)
!��������������� ������� ��� find_d_direction
use gen_mod
use lines_mod
external f_dwdz
real(8) f_d_dir,g,tt
complex(8) f_dwdz,z0
tt=-zarg(f_dwdz(z0+lines_d_dir_dr*cdexp(ii*g)))
if (tt<d0) tt=tt+pi2
f_d_dir=g-tt
end

subroutine get_lines_ds(dwdir0,dwzdir0,ds,z,kdir,dwdir,dwzdir,z12,z1,dzdzeta_z12,f_dwdzeta,f_dzdzeta)
!������ ����������� ���� ds ��� �������������� ��������� ����� ����
!dwdir0 - ��������� ������ ����������� ������� �������� � ������� ����� z (� ��������� zeta)
!dwzdir0 - ��������� ������ ����������� ������� �������� � ������� ����� z (� ��������� z)
!ds - ������� � �������� �������� ���� �������������� 
!z - ������� ����� � ��������� zeta
!kdir - 1 - �� ������, -1 - ������ ������
!dwdir - ��������� ������ ����������� ������� �������� � ����� z1 � ��������� zeta (��������)
!dwzdir - ��������� ������ ����������� ������� �������� � ����� z1 � ��������� z (��������)
!z12 - ������������� ����� �������������� � ��������� zeta (��������)
!z1 - ��������� ����� �������������� � ��������� zeta (��������)
!f_dwdzeta - ������� ������� dw/dzeta
!f_dzdzeta - ������� ������� dz/dzeta
use gen_mod
use lines_mod
external f_dwdzeta,f_dzdzeta
complex(8) f_dwdzeta,f_dzdzeta
complex(8) dwdir0,dwzdir0,z,dwdir,dwzdir,z12,dw,z1,dwz,dzdzeta_z12
real(8) ds,err_dir,err_dirz,adw,adwz
integer(4) kdir,j
j=0
do while(.true.)
    if (ds<lines_dsmin) then
        ds=lines_dsmin
        j=-9
    endif
    if (ds>lines_dsmax) then
        ds=lines_dsmax
        j=-9
    endif
    z12=z+dwdir0*ds*d5
	dw=f_dwdzeta(z12)
	dzdzeta_z12=f_dzdzeta(z12)
	dwz=dconjg(dw/dzdzeta_z12)*kdir
	dw=dconjg(dw)*kdir
	adw=cdabs(dw)
	adwz=cdabs(dwz)
	if (adw<1.0d-8.or.adwz<1.0d-8) then
		dwdir=c0  !������ � ����������� �����
		dwzdir=c0
		z1=z12
		exit
    else
		dwdir=dw/adw
		dwzdir=dwz/adwz
	endif
    err_dir=cdabs(dwdir0-dwdir)
	err_dirz=cdabs(dwzdir0-dwzdir)
    if ((err_dir<=lines_eps_dir).or.(err_dirz<=lines_eps_dir).or.(j==-9)) then
      z1=z+dwdir*ds
	  dw=f_dwdzeta(z1)
	  dwz=dconjg(dw/f_dzdzeta(z1))*kdir
	  dw=dconjg(dw)*kdir
      adw=cdabs(dw)
	  adwz=cdabs(dwz)
	  if (adw<1.0d-8.or.adwz<1.0d-8) then
	      dwdir=c0   !������ � ����������� �����
		  dwzdir=c0
		  exit
      else
		  dwdir=dw/adw
		  dwzdir=dwz/adwz
	  endif
      if (j==-9) exit
      err_dir=cdabs(dwdir0-dwdir)
	  err_dirz=cdabs(dwzdir0-dwzdir)
    endif  
    if(err_dir>lines_eps_dir.or.err_dirz>lines_eps_dir) then
        ds=ds/lines_kk
        if (j==1) then    
            j=-9
        else            
            j=-1
        endif
        cycle
    endif
    if(err_dir<lines_eps_dir*d5.and.err_dirz<lines_eps_dir*d5) then
        ds=ds*lines_kk
        if (j==-1) then  
            j=-9
        else            
            j=1
        endif
        cycle
    endif
    exit
enddo
end

subroutine find_line(z0,dir0,kdir,zl,k,nmax,f_dwdzeta,f_test_stop)
!���������� ����� ���� � ��������� zeta
!z0 - ��������� ����� � ��������� zeta
!dir0 - ��������� ����������� (� ��������� zeta)
!kdir - 1 - �� ������, -1 - ������ ������
!zl - �������� ������ ����� � ��������� zeta
!k - ��������� ������ � �������� ����� (����� ����� - 1)
!nmax - ������������ ������ � �������� �����
!f_dwdzeta - ������� ������� dw/dzeta
!f_test_stop - ������� ������� ��� �������� ��� ������� �� ��������� [logacal function lines_test_stop_empty(complex(8))]
             !������� �������� - ������� ����� �������������� � ���������� zeta
             !��������� .true. - �����������, .false. - �� ���������������
             !���� ��������������� �������� �� ��������� ���, �� ���������� � �������� ��������� lines_test_stop_empty
integer(4) kdir,k,nmax
complex(8) z0,zl(0:nmax)
complex(8), allocatable :: zlz(:)
real(8) dir0
external f_dwdzeta,f_dzdzeta_empty,f_test_stop
complex(8) f_dwdzeta,f_dzdzeta_empty
logical f_test_stop
allocate(zlz(nmax+1))
call find_line2(z0,z0,dir0,dir0,kdir,zl,zlz,k,nmax,f_dwdzeta,f_dzdzeta_empty,f_test_stop)
deallocate(zlz)
end

function f_dzdzeta_empty(z)
use gen_mod
complex(8) f_dzdzeta_empty,z, z1
z1 = z
f_dzdzeta_empty=c1
end

function lines_test_stop_empty(z,zz)
use gen_mod
complex(8) z, z1,zz
logical lines_test_stop_empty
z1 = z
z1 = zz
lines_test_stop_empty=.false.
end

subroutine find_line2(z0,zz0,dir0,dirz0,kdir,zl,zlz,k,nmax,f_dwdzeta,f_dzdzeta,f_test_stop)
!���������� ����� ���� � ��������� zeta � z
!z0 - ��������� ����� � ��������� zeta
!zz0 - ��������� ����� � ��������� z
!dir0 - ��������� ����������� (� ��������� zeta)
!dirz0 - ��������� ����������� (� ��������� z)
!kdir - 1 - �� ������, -1 - ������ ������
!zl - �������� ������ ����� � ��������� zeta
!zlz - �������� ������ ����� � ��������� z
!k - ��������� ������ � �������� ����� (����� ����� - 1)
!nmax - ������������ ������ � �������� �����
!f_dwdzeta - ������� ������� dw/dzeta
!f_dzdzeta - ������� ������� dz/dzeta
!f_test_stop - ������� ������� ��� �������� ��� ������� �� ��������� [logacal function lines_test_stop_empty(complex(8),complex(8))]
             !������� �������� - ������� ����� �������������� � ���������� zeta � z
             !��������� .true. - �����������, .false. - �� ���������������
             !���� ��������������� �������� �� ��������� ���, �� ���������� � �������� ��������� lines_test_stop_empty
use gen_mod
use lines_mod
integer(4) kdir,k,nmax
complex(8) z0,zl(0:nmax),zlz(0:nmax),z12,dwdir,dwdir0,darg1,darg2,zz0,dwzdir,dwzdir0,dzdzeta_z12
real(8) dir0,ds,ds_end,dsqrt2,adwdir,dirz0 
external f_dwdzeta,f_dzdzeta,f_test_stop
complex(8) f_dwdzeta,f_dzdzeta
logical f_test_stop !,qqq
dsqrt2=dsqrt(d2)
ds=lines_dsmax
zl=c0
zl(0)=z0
zlz(0)=zz0
call get_lines_ds(cdexp(ii*dir0),cdexp(ii*dirz0),ds,zl(0),kdir,dwdir,dwzdir,z12,zl(1),dzdzeta_z12,f_dwdzeta,f_dzdzeta)
!qqq=f_test_stop(zl(1))
ds_end=ds*d2
zlz(1)=zlz(0)+dzdzeta_z12*(zl(1)-zl(0))
k=1
do while(.true.)
    dwdir0=dwdir
	dwzdir0=dwzdir
    call get_lines_ds(dwdir0,dwzdir0,ds,zl(k),kdir,dwdir,dwzdir,z12,zl(k+1),dzdzeta_z12,f_dwdzeta,f_dzdzeta)
    zlz(k+1)=zlz(k)+dzdzeta_z12*(zl(k+1)-zl(k))
    k=k+1
	if (dwdir==d0) exit
	adwdir=cdabs(dwdir-dwdir0)
	if (adwdir>dsqrt2) exit
    if (k==nmax) exit
    if (cdabs(zl(k))>lines_rmax) exit
    if (k>5) then                    !��� ��������� �����
        darg1=zl(k)-z0
        darg1=darg1/cdabs(darg1)
        darg2=zl(k-1)-z0
        darg2=darg2/cdabs(darg2)
        if (cdabs(darg1-darg2)>dsqrt2) then !��������, ��� ���������� ����� ��������� ����� 
            zl(k)=z0
            zlz(k)=zlz(0)
            exit
        endif
        if (cdabs(zl(k)-z0)<ds_end) then !��������, ��� �������� � ��������� �����
            k=k+1
            zl(k)=z0
            zlz(k)=zlz(0)
            exit
        endif   
    endif
	if (f_test_stop(zl(k),zlz(k))) exit
enddo
end

subroutine rebuild_line(n,s,f,dsmin,dsmax,df,nmax,nf)
!������������� �������� ������� �� �������� ������������
use gen_mod
integer(4), parameter :: nfmax=2
integer(4) n         ! ����� ����� � �������
integer(4) nmax      !������������ ����� ����� � ����� �������
integer(4) nf        !����� �������
real(8) s(0:nmax)    !������� ��������
real(8) f(0:nmax,nf) !������� �������� �������
real(8) dsmin        !����������� ds
real(8) dsmax        !������������ ds
real(8) df(nf)       !����������� �� ������ ��������� ������� � �������� ������
integer(4) kdf(nfmax),i,i0,i1,k,n1,j
real(8) s2,f2,ds_prev,kk
real(8) s1(0:nmax),f1(0:nmax,nfmax)
real(8) bb_s(0:nmax,nfmax),cc_s(4,0:nmax,nfmax),bb_f(0:nmax,nfmax),cc_f(4,0:nmax,nfmax)
logical isend,dfcorrect,i1_finded

if (nf>nfmax) then
  write (*,*) "!!! nfmax error (subroutine rebuild_line)"
  stop
endif
kk=1.5d0 ! ����������� �� ����������� ���������� ds ����� ��������� ��������� 
!kk=d0
i0=0
i1=0
s1=d0
f1=d0
k=0
s1(k)=s(0)
forall(j=1:nf) f1(k,j)=f(0,j)
do while (i1<n)
  !���� ������� ������������
  kdf=0
  i0=i1
  do i=i0+1,n
    i1=i
	i1_finded=.false.
	do j=1,nf
      if (kdf(j)==0) then
	    if (f(i,j)>f(i-1,j)) then
	      kdf(j)=1
	    elseif (f(i,j)<f(i-1,j)) then
	      kdf(j)=-1
	    endif
	  elseif (kdf(j)>0) then
        if (f(i,j)<f(i-1,j)) then
	      i1=i-1
		  i1_finded=.true.
	      exit
	    endif
	  elseif (kdf(j)<0) then
        if (f(i,j)>f(i-1,j)) then
	      i1=i-1
		  i1_finded=.true.
	      exit
	    endif
	  endif
	enddo
	if (i1_finded) exit
  enddo
  !�������������
  n1=i1-i0
  do j=1,nf
    !call dcsakm(n1+1,s(i0:i1),f(i0:i1,j),bb_f(:,j),cc_f(:,:,j))
    !call dcsakm(n1+1,f(i0:i1,j),s(i0:i1),bb_s(:,j),cc_s(:,:,j))
    call spline_linear(n1+1,s(i0:i1),f(i0:i1,j),bb_f(:,j),cc_f(:,:,j))
    call spline_linear(n1+1,f(i0:i1,j),s(i0:i1),bb_s(:,j),cc_s(:,:,j))
  enddo
  isend=.false.
  do while (.not.isend)
    s2=s1(k)+dsmax
	if (s2>s(i1)) then
	  s2=s(i1)
	  isend=.true.
	endif
	!�������� ����������� df
	dfcorrect=.false.
	do j=1,nf
	  if (kdf(j)==0) cycle
      f2=dcsval(s2,n1,bb_f(:,j),cc_f(:,:,j))
	  if (dabs(f2-f1(k,j))>df(j)) then
	    f2=f1(k,j)+kdf(j)*df(j)
	    s2=dcsval(f2,n1,bb_s(:,j),cc_s(:,:,j))
		dfcorrect=.true.
      endif
	enddo
	if (dfcorrect) then
	  isend=.false.
	  !�������� ����������� kk
	  if ((k>0).and.(kk>d0)) then
	    ds_prev=(s1(k)-s1(k-1))*kk
	    if (s2-s1(k)>ds_prev) s2=s1(k)+ds_prev
      endif
	  !�������� �� ����������� ds
	  if (s2-s1(k)<dsmin) s2=s1(k)+dsmin
	  !�������� �� ���������� ����� ������� ������������
	  if (s2>s(i1)) then
	    s2=s(i1)
		isend=.true.
      endif
    endif
	!���������� �������� ������� �� ���������� s2
	if (.not.((s2-s1(k)<dsmin).and.isend)) k=k+1
	s1(k)=s2
    do j=1,nf
	  if (isend) then
	    f1(k,j)=f(i1,j)
	  else
	    f1(k,j)=dcsval(s2,n1,bb_f(:,j),cc_f(:,:,j))
	  endif
    enddo
	if ((k==nmax).and.(.not.isend)) then
	  isend=.true.
	  write(*,*)"!!! nmax error (subroutine rebuild_line)"
	endif
  enddo
enddo
n=k
s=s1
forall(j=1:nf) f(:,j)=f1(:,j)
end