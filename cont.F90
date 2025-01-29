! gfortran -o cont cont.F90 -ffree-line-length-none -I/home/xingwei/fftw/include/ -L/home/xingwei/fftw/lib/ -lfftw3

#define flow
!#define field

program main
implicit none

! constants
double complex one
parameter (one=(0.d0,1.d0))
double precision pi
parameter (pi=acos(-1.d0))

! dimensions of arrays, nc is the first dimension in spectral space
integer n1, n2, n3
parameter (n1=128, n2=128, n3=128)
integer nc
parameter (nc=n1/2+1)
integer n1o, n2o, n3o
parameter (n1o=4, n2o=4, n3o=4)
integer nco
parameter (nco=n1o/2+1)

! grids in physical space
double precision x1(n1), x2(n2), x3(n3)

! wavenumbers in spectral space
double precision kx(nc), ky(n2), kz(n3)

! velocity variables
double complex    u1(nc,n2,n3),  u2(nc,n2,n3),  u3(nc,n2,n3) ! velocity in spectral space
double precision  v1(n1,n2,n3),  v2(n1,n2,n3),  v3(n1,n2,n3) ! velocity in physical space
double complex u1o(nco,n2o,n3o), u2o(nco,n2o,n3o), u3o(nco,n2o,n3o) ! velocity of the low resolution

! field variables
double complex    a1(nc,n2,n3),  a2(nc,n2,n3),  a3(nc,n2,n3) ! field in spectral space
double precision  b1(n1,n2,n3),  b2(n1,n2,n3),  b3(n1,n2,n3) ! field in physical space
double complex a1o(nco,n2o,n3o), a2o(nco,n2o,n3o), a3o(nco,n2o,n3o) ! field of the low resolution

! for read in
double precision i1, i2, i3
double precision w1r, w1i, w2r, w2i, w3r, w3i

! for fftw
integer*8 pf, pb
double precision phys(n1,n2,n3)
double complex   spec(nc,n2,n3)
integer fftw_measure
parameter (fftw_measure=0)

integer j1, j2, j3, k1, k2, k3

call dfftw_plan_dft_r2c_3d(pf,n1,n2,n3,phys,spec,fftw_measure)
call dfftw_plan_dft_c2r_3d(pb,n1,n2,n3,spec,phys,fftw_measure)

! grids in physical space
do j1=1,n1
 x1(j1)=2*pi/n1*(j1-1)
enddo
do j2=1,n2
 x2(j2)=2*pi/n2*(j2-1)
enddo
do j3=1,n3
 x3(j3)=2*pi/n3*(j3-1)
enddo

! wavenumbers in spectral space
do k1=1,nc
 kx(k1)=k1-1
enddo
do k2=1,n2
 ky(k2)=k2-1
 if((k2-1).gt.n2/2) ky(k2)=k2-1-n2
enddo
do k3=1,n3
 kz(k3)=k3-1
 if((k3-1).gt.n3/2) kz(k3)=k3-1-n3
enddo

#ifdef flow
open(1,file='u',form='formatted')
#endif
#ifdef field
open(2,file='a',form='formatted')
#endif
do k3=1,n3o
 do k2=1,n2o
  do k1=1,nco
#ifdef flow
   read(1, '(3F6.1,6E15.6)') i1, i2, i3, w1r, w1i, w2r, w2i, w3r, w3i
   u1(k1,k2,k3)=w1r+one*w1i
   u2(k1,k2,k3)=w2r+one*w2i
   u3(k1,k2,k3)=w3r+one*w3i
#endif
#ifdef field
   read(2, '(3F6.1,6E15.6)') i1, i2, i3, w1r, w1i, w2r, w2i, w3r, w3i
   a1(k1,k2,k3)=w1r+one*w1i
   a2(k1,k2,k3)=w2r+one*w2i
   a3(k1,k2,k3)=w3r+one*w3i
#endif
  enddo
 enddo
enddo
#ifdef flow
close(1)
#endif
#ifdef field
close(2)
#endif

do k3=n3o+1,n3
 do k2=n2o+1,n2
  do k1=nco+1,nc
#ifdef flow
   u1(k1,k2,k3)=(0.d0,0.d0)
   u2(k1,k2,k3)=(0.d0,0.d0)
   u3(k1,k2,k3)=(0.d0,0.d0)
#endif
#ifdef field
   a1(k1,k2,k3)=(0.d0,0.d0)
   a2(k1,k2,k3)=(0.d0,0.d0)
   a3(k1,k2,k3)=(0.d0,0.d0)
#endif
  enddo
 enddo
enddo

#ifdef flow
spec=u1
call dfftw_execute_dft_c2r(pb,spec,phys)
v1=phys
spec=u2
call dfftw_execute_dft_c2r(pb,spec,phys)
v2=phys
spec=u3
call dfftw_execute_dft_c2r(pb,spec,phys)
v3=phys
#endif
#ifdef field
spec=a1
call dfftw_execute_dft_c2r(pb,spec,phys)
b1=phys
spec=a2
call dfftw_execute_dft_c2r(pb,spec,phys)
b2=phys
spec=a3
call dfftw_execute_dft_c2r(pb,spec,phys)
b3=phys
#endif

! output v in x-z plane
open(1,file='xz',form='formatted')
do j3=1,n3
 do j1=1,n1 
  write(1,'(6E15.6)') x1(j1), x3(j3), v1(j1,n2/2+1,j3), v2(j1,n2/2+1,j3), v3(j1,n2/2+1,j3), v1(j1,n2/2+1,j3)**2+v2(j1,n2/2+1,j3)**2+v3(j1,n2/2+1,j3)**2
 enddo
enddo
close(1)

end program main
