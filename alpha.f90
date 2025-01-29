program main
implicit none

double complex one
parameter (one=(0.d0,1.d0))
double precision pi
parameter (pi=acos(-1.d0))

integer n1, n2, n3
parameter (n1=128, n2=128, n3=128)
integer nc
parameter (nc=n1/2+1)

double precision  x1(n1), x2(n2), x3(n3)
double complex    a1(nc,n2,n3),  a2(nc,n2,n3),  a3(nc,n2,n3) ! field in spectral space
double precision  b1(n1,n2,n3),  b2(n1,n2,n3),  b3(n1,n2,n3) ! field in physical space
double precision  h1(n1,n2,n3),  h2(n1,n2,n3),  h3(n1,n2,n3)
double complex   a01(nc,n2,n3), a02(nc,n2,n3), a03(nc,n2,n3)
double precision b01(n1,n2,n3), b02(n1,n2,n3), b03(n1,n2,n3)
double precision  v1(n1,n2,n3),  v2(n1,n2,n3),  v3(n1,n2,n3)

integer*8 pf, pb
double precision phys(n1,n2,n3)
double complex   spec(nc,n2,n3)
integer fftw_measure
parameter (fftw_measure=0)

double precision rm, wnum

integer k1, k2, k3, j1, j2, j3
double precision i1, i2, i3, w1r, w1i, w2r, w2i, w3r, w3i
double precision y1, y2, y3

call dfftw_plan_dft_r2c_3d(pf,n1,n2,n3,phys,spec,fftw_measure)
call dfftw_plan_dft_c2r_3d(pb,n1,n2,n3,spec,phys,fftw_measure)

do j1=1,n1
 x1(j1)=2*pi/n1*(j1-1)
enddo
do j2=1,n2
 x2(j2)=2*pi/n2*(j2-1)
enddo
do j3=1,n3
 x3(j3)=2*pi/n3*(j3-1)
enddo

open(1,file='a',form='formatted')
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   read(1, '(3F6.1,6F15.6)') i1, i2, i3, w1r, w1i, w2r, w2i, w3r, w3i
   a1(k1,k2,k3)=w1r+one*w1i
   a2(k1,k2,k3)=w2r+one*w2i
   a3(k1,k2,k3)=w3r+one*w3i
  enddo
 enddo
enddo
close(1)
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   a01(k1,k2,k3)=(0.d0,0.d0)
   a02(k1,k2,k3)=(0.d0,0.d0)
   a03(k1,k2,k3)=(0.d0,0.d0)
  enddo
 enddo
enddo
do k3=1,2
 do k2=1,2
  do k1=1,2
   a01(k1,k2,k3)=a1(k1,k2,k3)
   a02(k1,k2,k3)=a2(k1,k2,k3)
   a03(k1,k2,k3)=a3(k1,k2,k3)
  enddo
 enddo
enddo
spec=a01
call dfftw_execute_dft_c2r(pb,spec,phys)
b01=phys
spec=a02
call dfftw_execute_dft_c2r(pb,spec,phys)
b02=phys
spec=a03
call dfftw_execute_dft_c2r(pb,spec,phys)
b03=phys

open(1,file='b',form='formatted')
open(2,file='h',form='formatted')
do j3=1,n3
 do j2=1,n2
  do j1=1,n1
   read(1,'(6F15.6)'), y1, y2, y3, b1(j1,j2,j3), b2(j1,j2,j3), b3(j1,j2,j3)
   read(2,'(6F15.6)'), y1, y2, y3, h1(j1,j2,j3), h2(j1,j2,j3), h3(j1,j2,j3)
  enddo
 enddo
enddo
close(1)
close(2)

do j3=1,n3
 do j2=1,n2
  do j1=1,n1
   v1(j1,j2,j3)=-sin(wnum*x2(j2))+cos(wnum*x3(j3))
   v2(j1,j2,j3)= cos(wnum*x1(j1))-sin(wnum*x3(j3))
   v3(j1,j2,j3)=-sin(wnum*x1(j1))+cos(wnum*x2(j2))
  enddo
 enddo
enddo

do j3=1,n3
 do j2=1,n2
  do j1=1,n1
   h1(j1,j2,j3)=v2(j1,j2,j3)*(b3(j1,j2,j3)-b03(j1,j2,j3))-v3(j1,j2,j3)*(b2(j1,j2,j3)-b02(j1,j2,j3))
   h2(j1,j2,j3)=v3(j1,j2,j3)*(b1(j1,j2,j3)-b01(j1,j2,j3))-v1(j1,j2,j3)*(b3(j1,j2,j3)-b03(j1,j2,j3))
   h3(j1,j2,j3)=v1(j1,j2,j3)*(b2(j1,j2,j3)-b02(j1,j2,j3))-v2(j1,j2,j3)*(b1(j1,j2,j3)-b01(j1,j2,j3))
   write(1,'(9E15.6)') b01(j1,j2,j3), h1(j1,j2,j3), h1(j1,j2,j3)/b01(j1,j2,j3), &
   &                   b02(j1,j2,j3), h2(j1,j2,j3), h2(j1,j2,j3)/b02(j1,j2,j3), &
   &                   b03(j1,j2,j3), h3(j1,j2,j3), h3(j1,j2,j3)/b03(j1,j2,j3)
  enddo
 enddo
enddo
close(1)

rm=10.d0
wnum=40.d0
write(6,*) -rm/wnum
stop
end program main
