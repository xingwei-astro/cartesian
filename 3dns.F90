! gfortran -fbounds-check -ffree-line-length-none -o 3dns 3dns.F90 -I$HOME/fftw/include/ -L$HOME/fftw/lib/ -lfftw3 -lm
! 3D rotating MHD problem, decaying and forced
! dv/dt = -grad p + (v x cv) + 1/Re*(del^2 v) + 2/Ro*(v x e_3) + Va^2*(cb x (B0 + b))
! db/dt = curl(v x (B0 + b)) + 1/Rm*(del^2 b)
! length normalized with box size l, velocity with initial mean velocity u0=sqrt(u^2), time with l/u0, field with B0
! Re=u0*l/nu, Ro=u0/(Omega*l), Va=B0/sqrt(rho*mu)/u0, Rm=u0*l/eta

! rotating flow driven by periodic helical forcing with amplitude af, wavenumber kf and frequency freq
! dv/dt = -grad p + (v x cv) + Ek*(del^2 v) + 2v x e_3 + f
! length normalized with box size l, time with rotation rate Omega^-1, and velocity with Omega*l

#define flow
#define partc
#define field
!#define kinematic  !this flag cannot co-exist with the flag #define flow

module global
implicit none
! constants
double complex one
parameter (one=(0.d0,1.d0))
double precision pi
parameter (pi=acos(-1.d0))

! dimensionless parameters
double precision Re, Ro, Va, Rm, Ek
parameter (Re=1.d3, Ro=1.d-1, Va=1.d-1, Rm=1.d3, Ek=1.d-3)

! driving force
double precision af, af1, af2, af3
integer kf1, kf2, kf3
double precision kf
double precision freq
double precision argf1, argf2, argf3
double precision f1r, f1i, f2r, f2i, f3r, f3i
parameter (af=1.d-3, kf1=1, kf2=1, kf3=1, freq=-2.d0/sqrt(3.d0))

! imposed field
double precision B01, B02, B03
double precision alpha !angle between rotational and magnetic axes
parameter (alpha=45.d0/180.d0*pi)
parameter (B01=sin(alpha), B02=0.d0, B03=cos(alpha))

! dimensions of arrays, nc is the first dimension in spectral space
integer n1, n2, n3
parameter (n1=128, n2=128, n3=128)
integer nc
parameter (nc=n1/2+1)

! time step
integer nt, nstep
parameter (nt=10, nstep=10)
double precision dt, t
parameter (dt=1.d-2)

! grids in physical space
double precision x1(n1), x2(n2), x3(n3)

! wavenumbers in spectral space
double precision kx(nc), ky(n2), kz(n3)

! velocity variables
double complex    u1(nc,n2,n3),  u2(nc,n2,n3),  u3(nc,n2,n3) ! velocity in spectral space
double precision  v1(n1,n2,n3),  v2(n1,n2,n3),  v3(n1,n2,n3) ! velocity in physical space
double complex   cu1(nc,n2,n3), cu2(nc,n2,n3), cu3(nc,n2,n3) ! vorticity in spectral space
double precision cv1(n1,n2,n3), cv2(n1,n2,n3), cv3(n1,n2,n3) ! vorticity in physical space

! field variables
double complex    a1(nc,n2,n3),  a2(nc,n2,n3),  a3(nc,n2,n3) ! field in spectral space
double precision  b1(n1,n2,n3),  b2(n1,n2,n3),  b3(n1,n2,n3) ! field in physical space
double complex   ca1(nc,n2,n3), ca2(nc,n2,n3), ca3(nc,n2,n3) ! current in spectral space
double precision cb1(n1,n2,n3), cb2(n1,n2,n3), cb3(n1,n2,n3) ! current in physical space

! for fftw
integer*8 pf, pb
double precision phys(n1,n2,n3)
double complex   spec(nc,n2,n3)
integer fftw_measure
parameter (fftw_measure=0)
end module global

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global
implicit none
integer j1, j2, j3, k1, k2, k3
integer it
double precision t0
double complex u1o(nc,n2,n3), u2o(nc,n2,n3), u3o(nc,n2,n3) ! velocity at last time layer
double complex a1o(nc,n2,n3), a2o(nc,n2,n3), a3o(nc,n2,n3) ! field at last time layer
double complex e1(nc,n2,n3), e2(nc,n2,n3), e3(nc,n2,n3) ! force in spectral space
double complex g1(nc,n2,n3), g2(nc,n2,n3), g3(nc,n2,n3) ! induction in spectral space
double complex e1_1(nc,n2,n3), e2_1(nc,n2,n3), e3_1(nc,n2,n3)
double complex e1_2(nc,n2,n3), e2_2(nc,n2,n3), e3_2(nc,n2,n3)
double complex g1_1(nc,n2,n3), g2_1(nc,n2,n3), g3_1(nc,n2,n3)
double complex g1_2(nc,n2,n3), g2_2(nc,n2,n3), g3_2(nc,n2,n3)
double precision enrv1, enrv2, enrv3, enrv, enrw1, enrw2, enrw3, enrw, divv 
double precision enrb1, enrb2, enrb3, enrb, enrj1, enrj2, enrj3, enrj, divb
double precision enrbo1, enrbo2, enrbo3, enrbo
integer wnum1, wnum2, wnum3
double precision freq1, freq2, freq3
double precision px1_0, px2_0, px3_0, & ! 0 center, 1 and 2 in x, 3 and 4 in y, 5 and 6 in z
&                px1_1, px2_1, px3_1, px1_2, px2_2, px3_2, & 
&                px1_3, px2_3, px3_3, px1_4, px2_4, px3_4, &
&                px1_5, px2_5, px3_5, px1_6, px2_6, px3_6

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

! driving force 
kf=sqrt(dble(kf1)**2+dble(kf2)**2+dble(kf3)**2)
af1=af*sqrt(dble(kf2)**2+dble(kf3)**2)/sqrt(2.d0)/kf
af2=af*sqrt(dble(kf3)**2+dble(kf1)**2)/sqrt(2.d0)/kf
af3=af*sqrt(dble(kf1)**2+dble(kf2)**2)/sqrt(2.d0)/kf
argf1=0.d0
argf2=argf1+pi-acos(dble(kf1)*dble(kf2)/sqrt((dble(kf2)**2+dble(kf3)**2)*(dble(kf3)**2+dble(kf1)**2)))
argf3=argf1+pi+acos(dble(kf3)*dble(kf1)/sqrt((dble(kf1)**2+dble(kf2)**2)*(dble(kf2)**2+dble(kf3)**2)))
f1r=af1*cos(argf1)
f1i=af1*sin(argf1)
f2r=af2*cos(argf2)
f2i=af2*sin(argf2)
f3r=af3*cos(argf3)
f3i=af3*sin(argf3)

! parameters of helical wave
#ifdef kinematic
wnum1=1
wnum2=1
wnum3=1
freq1=2*pi/1.d2
freq2=2*pi/1.d2
freq3=2*pi/1.d2
#endif

! particle initial position
#ifdef partc
px1_0=pi
px2_0=pi
px3_0=pi
px1_1=pi+1.d-2
px2_1=pi
px3_1=pi
px1_2=pi-1.d-2
px2_2=pi
px3_2=pi
px1_3=pi
px2_3=pi+1.d-2
px3_3=pi
px1_4=pi
px2_4=pi-1.d-2
px3_4=pi
px1_5=pi
px2_5=pi
px3_5=pi+1.d-2
px1_6=pi
px2_6=pi
px3_6=pi-1.d-2
#endif

! initialise
call init(t0)

! open files for monitor
#ifdef flow
open(1,file='flow',form='formatted')
#endif
#ifdef field
open(2,file='field',form='formatted')
#endif
#ifdef partc
open(3,file='partc',form='formatted')
#endif

! time stepping
do it=1,nt
 t=t0+it*dt

! given flow of helical wave
#ifdef kinematic
 do j3=1,n3
  do j2=1,n2
   do j1=1,n1
! one wave
    v1(j1,j2,j3)=0
    v2(j1,j2,j3)=cos(wnum1*x1(j1)-freq1*t)
    v3(j1,j2,j3)=-sin(wnum1*x1(j1)-freq1*t)
! two waves
    v1(j1,j2,j3)=v1(j1,j2,j3)-sin(wnum2*x2(j2)-freq2*t)
    v3(j1,j2,j3)=v3(j1,j2,j3)+cos(wnum2*x2(j2)-freq2*t)
! three waves
!    v1(j1,j2,j3)=v1(j1,j2,j3)+cos(wnum3*x3(j3)-freq3*t)
!    v2(j1,j2,j3)=v2(j1,j2,j3)-sin(wnum3*x3(j3)-freq3*t)
   enddo
  enddo
 enddo
#endif

! save velocity and field of the last time layer
#ifdef flow
 u1o=u1
 u2o=u2
 u3o=u3
#endif
#ifdef field
 a1o=a1
 a2o=a2
 a3o=a3
#endif

! for kinematic dynamo
#ifdef kinematic
 call energy(a1,a2,a3,enrbo1,enrbo2,enrbo3,enrbo)
#endif

! first call nonlin to calculate forces e_1 and g_1
 call nonlin(e1_1,e2_1,e3_1,g1_1,g2_1,g3_1)

! first call step to update u and a with e_1 and g_1
 call step(u1o,u2o,u3o,a1o,a2o,a3o,e1_1,e2_1,e3_1,g1_1,g2_1,g3_1)

! second call nonlin to calculate force e_2 and g_2
 call nonlin(e1_2,e2_2,e3_2,g1_2,g2_2,g3_2)

#ifdef flow
 e1=(e1_1+e1_2)/2
 e2=(e2_1+e2_2)/2
 e3=(e3_1+e3_2)/2
#endif
#ifdef field
 g1=(g1_1+g1_2)/2
 g2=(g2_1+g2_2)/2
 g3=(g3_1+g3_2)/2
#endif

! second call step to update u and a with e and g
 call step(u1o,u2o,u3o,a1o,a2o,a3o,e1,e2,e3,g1,g2,g3)

! calculate particle position
#ifdef partc
 call particle(px1_0,px2_0,px3_0)
 call particle(px1_1,px2_1,px3_1)
 call particle(px1_2,px2_2,px3_2)
 call particle(px1_3,px2_3,px3_3)
 call particle(px1_4,px2_4,px3_4)
 call particle(px1_5,px2_5,px3_5)
 call particle(px1_6,px2_6,px3_6)
#endif

! monitor
if(mod(it,nt/nstep).eq.0) then
#ifdef flow
 call energy(u1,u2,u3,enrv1,enrv2,enrv3,enrv)
 call curl(u1,u2,u3,cu1,cu2,cu3)
 call energy(cu1,cu2,cu3,enrw1,enrw2,enrw3,enrw)
 call diverg(u1,u2,u3,divv)
 write(1,'(10E15.6)') t, enrv1, enrv2, enrv3, enrv, enrw1, enrw2, enrw3, enrw, divv
#ifdef partc
 write(3,'(22F6.3)') t, px1_0, px2_0, px3_0, &
 &                      px1_1, px2_1, px3_1, px1_2, px2_2, px3_2, &
 &                      px1_3, px2_3, px3_3, px1_4, px2_4, px3_4, &
 &                      px1_5, px2_5, px3_5, px1_6, px2_6, px3_6
#endif
#endif
#ifdef field
 call energy(a1,a2,a3,enrb1,enrb2,enrb3,enrb)
 call curl(a1,a2,a3,ca1,ca2,ca3)
 call energy(ca1,ca2,ca3,enrj1,enrj2,enrj3,enrj)
 call diverg(a1,a2,a3,divb)
#ifdef kinematic
 write(2,'(3E15.6)') t, log(enrb/enrbo)/dt, divb
#else
 write(2,'(10E15.6)') t, enrb1, enrb2, enrb3, enrb, enrj1, enrj2, enrj3, enrj, divb
#endif
#endif
endif

! normalisation for kinematic dynamo
#ifdef kinematic
 a1=a1/sqrt(enrb)
 a2=a2/sqrt(enrb)
 a3=a3/sqrt(enrb)
#endif

! end of this time step
if(abs(u1(2,2,2)).gt.1.d50) stop
if(abs(u2(2,2,2)).gt.1.d50) stop
if(abs(u3(2,2,2)).gt.1.d50) stop
if(abs(a1(2,2,2)).gt.1.d50) stop
if(abs(a2(2,2,2)).gt.1.d50) stop
if(abs(a3(2,2,2)).gt.1.d50) stop
enddo
#ifdef flow
close(1)
#endif
#ifdef field
close(2)
#endif
#ifdef partc
close(3)
#endif

! output
call output

call dfftw_destroy_plan(pf)
call dfftw_destroy_plan(pb)
end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nonlin(e1,e2,e3,g1,g2,g3)
use global
implicit none
double complex e1(nc,n2,n3), e2(nc,n2,n3), e3(nc,n2,n3)
double complex g1(nc,n2,n3), g2(nc,n2,n3), g3(nc,n2,n3)
double precision f1(n1,n2,n3), f2(n1,n2,n3), f3(n1,n2,n3)
double precision h1(n1,n2,n3), h2(n1,n2,n3), h3(n1,n2,n3)
integer j1, j2, j3, k1, k2, k3
! calculate v and b in physical space
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
! calculate cu and ca in spectral space
#ifdef flow
call curl(u1,u2,u3,cu1,cu2,cu3)
#ifdef field
call curl(a1,a2,a3,ca1,ca2,ca3)
#endif
#endif
! calculate cv and cb in physical space
#ifdef flow
spec=cu1
call dfftw_execute_dft_c2r(pb,spec,phys)
cv1=phys
spec=cu2
call dfftw_execute_dft_c2r(pb,spec,phys)
cv2=phys
spec=cu3
call dfftw_execute_dft_c2r(pb,spec,phys)
cv3=phys
#ifdef field
spec=ca1
call dfftw_execute_dft_c2r(pb,spec,phys)
cb1=phys
spec=ca2
call dfftw_execute_dft_c2r(pb,spec,phys)
cb2=phys
spec=ca3
call dfftw_execute_dft_c2r(pb,spec,phys)
cb3=phys
#endif
#endif
! calculate f and h in physical space
do j3=1,n3
 do j2=1,n2
  do j1=1,n1
#ifdef flow
! inertial force (v x cv)
   f1(j1,j2,j3)=v2(j1,j2,j3)*cv3(j1,j2,j3)-v3(j1,j2,j3)*cv2(j1,j2,j3)
   f2(j1,j2,j3)=v3(j1,j2,j3)*cv1(j1,j2,j3)-v1(j1,j2,j3)*cv3(j1,j2,j3)
   f3(j1,j2,j3)=v1(j1,j2,j3)*cv2(j1,j2,j3)-v2(j1,j2,j3)*cv1(j1,j2,j3)
!   f1(j1,j2,j3)=0.d0
!   f2(j1,j2,j3)=0.d0
!   f3(j1,j2,j3)=0.d0
! Coriolis force 2*(v x e_3)
   f1(j1,j2,j3)=f1(j1,j2,j3)+2.d0*v2(j1,j2,j3)
   f2(j1,j2,j3)=f2(j1,j2,j3)-2.d0*v1(j1,j2,j3)
#ifdef field
! Lorentz force (cb x b) for dynamo
   f1(j1,j2,j3)=f1(j1,j2,j3)+(b3(j1,j2,j3)*cb2(j1,j2,j3)-b2(j1,j2,j3)*cb3(j1,j2,j3))
   f2(j1,j2,j3)=f2(j1,j2,j3)+(b1(j1,j2,j3)*cb3(j1,j2,j3)-b3(j1,j2,j3)*cb1(j1,j2,j3))
   f3(j1,j2,j3)=f3(j1,j2,j3)+(b2(j1,j2,j3)*cb1(j1,j2,j3)-b1(j1,j2,j3)*cb2(j1,j2,j3))
! Lorentz force Va^2*(cb x (B0 + b)) for imposed field
!   f1(j1,j2,j3)=f1(j1,j2,j3)+Va**2*((B03+b3(j1,j2,j3))*cb2(j1,j2,j3)-(B02+b2(j1,j2,j3))*cb3(j1,j2,j3))
!   f2(j1,j2,j3)=f2(j1,j2,j3)+Va**2*((B01+b1(j1,j2,j3))*cb3(j1,j2,j3)-(B03+b3(j1,j2,j3))*cb1(j1,j2,j3))
!   f3(j1,j2,j3)=f3(j1,j2,j3)+Va**2*((B02+b2(j1,j2,j3))*cb1(j1,j2,j3)-(B01+b1(j1,j2,j3))*cb2(j1,j2,j3))
#endif
! driving force
!   f1(j1,j2,j3)=f1(j1,j2,j3)+af*(cos(kf2*x2(j2)-freq*t)+sin(kf3*x3(j3)-freq*t))
!   f2(j1,j2,j3)=f2(j1,j2,j3)+af*(cos(kf3*x3(j3)-freq*t)+sin(kf1*x1(j1)-freq*t))
!   f3(j1,j2,j3)=f3(j1,j2,j3)+af*(cos(kf1*x1(j1)-freq*t)+sin(kf2*x2(j2)-freq*t))
   f1(j1,j2,j3)=f1(j1,j2,j3)+f1r*cos(kf1*x1(j1)+kf2*x2(j2)+kf3*x3(j3)-freq*t)-f1i*sin(kf1*x1(j1)+kf2*x2(j2)+kf3*x3(j3)-freq*t)
   f2(j1,j2,j3)=f2(j1,j2,j3)+f2r*cos(kf1*x1(j1)+kf2*x2(j2)+kf3*x3(j3)-freq*t)-f2i*sin(kf1*x1(j1)+kf2*x2(j2)+kf3*x3(j3)-freq*t)
   f3(j1,j2,j3)=f3(j1,j2,j3)+f3r*cos(kf1*x1(j1)+kf2*x2(j2)+kf3*x3(j3)-freq*t)-f3i*sin(kf1*x1(j1)+kf2*x2(j2)+kf3*x3(j3)-freq*t)
#endif
#ifdef field
! induction term (v x b) for dynamo
   h1(j1,j2,j3)=v2(j1,j2,j3)*b3(j1,j2,j3)-v3(j1,j2,j3)*b2(j1,j2,j3)
   h2(j1,j2,j3)=v3(j1,j2,j3)*b1(j1,j2,j3)-v1(j1,j2,j3)*b3(j1,j2,j3)
   h3(j1,j2,j3)=v1(j1,j2,j3)*b2(j1,j2,j3)-v2(j1,j2,j3)*b1(j1,j2,j3)
! induction term (v x (B0 + b)) for imposed field
!   h1(j1,j2,j3)=v2(j1,j2,j3)*(B03+b3(j1,j2,j3))-v3(j1,j2,j3)*(B02+b2(j1,j2,j3))
!   h2(j1,j2,j3)=v3(j1,j2,j3)*(B01+b1(j1,j2,j3))-v1(j1,j2,j3)*(B03+b3(j1,j2,j3))
!   h3(j1,j2,j3)=v1(j1,j2,j3)*(B02+b2(j1,j2,j3))-v2(j1,j2,j3)*(B01+b1(j1,j2,j3))
! test for alpha effect
!   h1(j1,j2,j3)=0.d0
!   h2(j1,j2,j3)=v3(j1,j2,j3)
!   h3(j1,j2,j3)=-v2(j1,j2,j3)
#endif
  enddo
 enddo
enddo
! calculate e and g in spectral space
#ifdef flow
phys=f1
call dfftw_execute_dft_r2c(pf,phys,spec)
e1=spec/n1/n2/n3
phys=f2
call dfftw_execute_dft_r2c(pf,phys,spec)
e2=spec/n1/n2/n3
phys=f3
call dfftw_execute_dft_r2c(pf,phys,spec)
e3=spec/n1/n2/n3
#endif
#ifdef field
phys=h1
call dfftw_execute_dft_r2c(pf,phys,spec)
g1=spec/n1/n2/n3
phys=h2
call dfftw_execute_dft_r2c(pf,phys,spec)
g2=spec/n1/n2/n3
phys=h3
call dfftw_execute_dft_r2c(pf,phys,spec)
g3=spec/n1/n2/n3
#endif
end subroutine nonlin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine step(u1o,u2o,u3o,a1o,a2o,a3o,e1,e2,e3,g1,g2,g3)
use global
implicit none
double complex u1o(nc,n2,n3), u2o(nc,n2,n3), u3o(nc,n2,n3)
double complex a1o(nc,n2,n3), a2o(nc,n2,n3), a3o(nc,n2,n3)
double complex e1(nc,n2,n3), e2(nc,n2,n3), e3(nc,n2,n3)
double complex g1(nc,n2,n3), g2(nc,n2,n3), g3(nc,n2,n3)
integer j1, j2, j3, k1, k2, k3
double precision kk
double complex w ! w=(k dot e)/k^2 to calculate (-grad p+f)
double complex cg1(nc,n2,n3), cg2(nc,n2,n3), cg3(nc,n2,n3) ! curl(vxb)
#ifdef field
call curl(g1,g2,g3,cg1,cg2,cg3)
#endif
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   kk=kx(k1)*kx(k1)+ky(k2)*ky(k2)+kz(k3)*kz(k3)
#ifdef flow
   if(k1.eq.1 .and. k2.eq.1 .and. k3.eq.1) then
    w=0
   else
    w=(kx(k1)*e1(k1,k2,k3)+ky(k2)*e2(k1,k2,k3)+kz(k3)*e3(k1,k2,k3))/kk
   endif
   e1(k1,k2,k3)=e1(k1,k2,k3)-kx(k1)*w
   e2(k1,k2,k3)=e2(k1,k2,k3)-ky(k2)*w
   e3(k1,k2,k3)=e3(k1,k2,k3)-kz(k3)*w
   u1(k1,k2,k3)=((1-dt*kk/2*Ek)*u1o(k1,k2,k3)+dt*e1(k1,k2,k3))/(1+dt*kk/2*Ek)
   u2(k1,k2,k3)=((1-dt*kk/2*Ek)*u2o(k1,k2,k3)+dt*e2(k1,k2,k3))/(1+dt*kk/2*Ek)
   u3(k1,k2,k3)=((1-dt*kk/2*Ek)*u3o(k1,k2,k3)+dt*e3(k1,k2,k3))/(1+dt*kk/2*Ek)
#endif
#ifdef field
   a1(k1,k2,k3)=((1-dt*kk/2/Rm)*a1o(k1,k2,k3)+dt*cg1(k1,k2,k3))/(1+dt*kk/2/Rm)
   a2(k1,k2,k3)=((1-dt*kk/2/Rm)*a2o(k1,k2,k3)+dt*cg2(k1,k2,k3))/(1+dt*kk/2/Rm)
   a3(k1,k2,k3)=((1-dt*kk/2/Rm)*a3o(k1,k2,k3)+dt*cg3(k1,k2,k3))/(1+dt*kk/2/Rm)
#endif
  enddo
 enddo
enddo
end subroutine step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init(t0)
use global
implicit none
integer ic ! ic=0 zero; ic=1 certain mode; ic=2 random; ic=3 Gaussian; ic=4 sin*sin; ic=5: readin
double precision t0
integer j1, j2, j3, k1, k2, k3
double precision x, sigma, enrv1, enrv2, enrv3, enrv
double precision i1, i2, i3
double precision w1r, w1i, w2r, w2i, w3r, w3i
ic=2

t0=0.d0
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
    u1(k1,k2,k3)=0
    u2(k1,k2,k3)=0
    u3(k1,k2,k3)=0
   cu1(k1,k2,k3)=0
   cu2(k1,k2,k3)=0
   cu3(k1,k2,k3)=0
    a1(k1,k2,k3)=0
    a2(k1,k2,k3)=0
    a3(k1,k2,k3)=0
   ca1(k1,k2,k3)=0
   ca2(k1,k2,k3)=0
   ca3(k1,k2,k3)=0
  enddo
 enddo
enddo
do j3=1,n3
 do j2=1,n2
  do j1=1,n1
    v1(j1,j2,j3)=0
    v2(j1,j2,j3)=0
    v3(j1,j2,j3)=0
   cv1(j1,j2,j3)=0
   cv2(j1,j2,j3)=0
   cv3(j1,j2,j3)=0
    b1(j1,j2,j3)=0
    b2(j1,j2,j3)=0
    b3(j1,j2,j3)=0
   cb1(j1,j2,j3)=0
   cb2(j1,j2,j3)=0
   cb3(j1,j2,j3)=0
  enddo
 enddo
enddo
a1(1,2,2)=1.d-3
a2(2,1,2)=1.d-3
a3(2,2,1)=1.d-3

if(ic.eq.1) then
t0=0.d0
#ifdef flow
k1=2
k2=2
k3=2
u1(k1,k2,k3)=1.d-3
u2(k1,k2,k3)=u1(k1,k2,k3)
u3(k1,k2,k3)=-(kx(k1)*u1(k1,k2,k3)+ky(k2)*u2(k1,k2,k3))/kz(k3)
#endif
#ifdef field
k1=2
k2=2
k3=2
a1(k1,k2,k3)=1.d-3
a2(k1,k2,k3)=a1(k1,k2,k3)
a3(k1,k2,k3)=-(kx(k1)*a1(k1,k2,k3)+ky(k2)*a2(k1,k2,k3))/kz(k3)
#endif
endif

if(ic.eq.2) then
t0=0.d0
call random_seed()
#ifdef flow
do j3=1,n3
 do j2=1,n2
  do j1=1,n1   
   call random_number(x)
   v1(j1,j2,j3)=1.d-3*x
   call random_number(x)
   v2(j1,j2,j3)=1.d-3*x
  enddo
 enddo
enddo
phys=v1
call dfftw_execute_dft_r2c(pf,phys,spec)
u1=spec/n1/n2/n3
phys=v2
call dfftw_execute_dft_r2c(pf,phys,spec)
u2=spec/n1/n2/n3
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   if(k3.eq.1) then
    u3(k1,k2,k3)=0
   else
    u3(k1,k2,k3)=-(kx(k1)*u1(k1,k2,k3)+ky(k2)*u2(k1,k2,k3))/kz(k3)
   endif
  enddo
 enddo
enddo
spec=u3
call dfftw_execute_dft_c2r(pb,spec,phys)
v3=phys
#endif
#ifdef field
do j3=1,n3
 do j2=1,n2
  do j1=1,n1
   call random_number(x)
   b1(j1,j2,j3)=1.d-6*x
   call random_number(x)
   b2(j1,j2,j3)=1.d-6*x
  enddo
 enddo
enddo
phys=b1
call dfftw_execute_dft_r2c(pf,phys,spec)
a1=spec/n1/n2/n3
phys=b2
call dfftw_execute_dft_r2c(pf,phys,spec)
a2=spec/n1/n2/n3
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   if(k3.eq.1) then
    a3(k1,k2,k3)=0
   else
    a3(k1,k2,k3)=-(kx(k1)*a1(k1,k2,k3)+ky(k2)*a2(k1,k2,k3))/kz(k3)
   endif
  enddo
 enddo
enddo
spec=a3
call dfftw_execute_dft_c2r(pb,spec,phys)
b3=phys
#endif
endif

if(ic.eq.3) then
t0=0.d0
#ifdef flow
sigma=pi
do j3=1,n3
 do j2=1,n2
  do j1=1,n1
   v1(j1,j2,j3)=exp(-((x2(j2)-pi)**2+(x3(j3)-pi)**2)/sigma**2)
   v2(j1,j2,j3)=exp(-((x3(j3)-pi)**2+(x1(j1)-pi)**2)/sigma**2)
   v3(j1,j2,j3)=exp(-((x1(j1)-pi)**2+(x2(j2)-pi)**2)/sigma**2)
  enddo
 enddo
enddo
phys=v1
call dfftw_execute_dft_r2c(pf,phys,spec)
u1=spec/n1/n2/n3
phys=v2
call dfftw_execute_dft_r2c(pf,phys,spec)
u2=spec/n1/n2/n3
phys=v3
call dfftw_execute_dft_r2c(pf,phys,spec)
u3=spec/n1/n2/n3
u1(1,1,1)=0
u2(1,1,1)=0
u3(1,1,1)=0
spec=u1
call dfftw_execute_dft_c2r(pb,spec,phys)
v1=phys
spec=u2
call dfftw_execute_dft_c2r(pb,spec,phys)
v2=phys
spec=u3
call dfftw_execute_dft_c2r(pb,spec,phys)
v3=phys
call energy(u1,u2,u3,enrv1,enrv2,enrv3,enrv)
u1=u1/sqrt(enrv)
u2=u2/sqrt(enrv)
u3=u3/sqrt(enrv)
#endif
endif

if(ic.eq.4) then
t0=0.d0
#ifdef flow
k1=1
k2=1
k3=1
do j3=1,n3
 do j2=1,n2
  do j1=1,n1
   v1(j1,j2,j3)=sin(k2*x2(j2))*sin(k3*x3(j3))
   v2(j1,j2,j3)=sin(k3*x3(j3))*sin(k1*x1(j1))
   v3(j1,j2,j3)=sin(k1*x1(j1))*sin(k2*x2(j2))
  enddo
 enddo
enddo
phys=v1
call dfftw_execute_dft_r2c(pf,phys,spec)
u1=spec/n1/n2/n3
phys=v2
call dfftw_execute_dft_r2c(pf,phys,spec)
u2=spec/n1/n2/n3
phys=v3
call dfftw_execute_dft_r2c(pf,phys,spec)
u3=spec/n1/n2/n3
call energy(u1,u2,u3,enrv1,enrv2,enrv3,enrv)
u1=u1/sqrt(enrv)
u2=u2/sqrt(enrv)
u3=u3/sqrt(enrv)
#endif
endif

if(ic.eq.5) then
open(1,file='time',form='formatted')
read(1,*) t0
close(1)
#ifdef flow
open(1,file='u',form='formatted')
#endif
#ifdef field
open(2,file='a',form='formatted')
#endif
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
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
endif
end subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output
use global
implicit none
integer j1, j2, j3, k1, k2, k3
double precision h1(n1,n2,n3), h2(n1,n2,n3), h3(n1,n2,n3)

open(1,file='time',form='formatted')
write(1,*) t
close(1)

#ifdef flow
open(1,file='u',form='formatted')
open(2,file='cu',form='formatted')
#endif
#ifdef field
open(3,file='a',form='formatted')
open(4,file='ca',form='formatted')
#endif
! calculate u in case of a given flow for kinematic dynamo
#ifdef kinematic
phys=v1
call dfftw_execute_dft_r2c(pf,phys,spec)
u1=spec/n1/n2/n3
phys=v2
call dfftw_execute_dft_r2c(pf,phys,spec)
u2=spec/n1/n2/n3
phys=v3
call dfftw_execute_dft_r2c(pf,phys,spec)
u3=spec/n1/n2/n3
#endif
! calculate cu and ca
#ifdef flow
call curl(u1,u2,u3,cu1,cu2,cu3)
#endif
#ifdef field
call curl(a1,a2,a3,ca1,ca2,ca3)
#endif
! output u, a, cu and ca
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
#ifdef flow
   write(1, '(3F6.1,6E15.6)') kx(k1), ky(k2), kz(k3), &
   &  dble(u1(k1,k2,k3)), imag(u1(k1,k2,k3)), &
   &  dble(u2(k1,k2,k3)), imag(u2(k1,k2,k3)), &
   &  dble(u3(k1,k2,k3)), imag(u3(k1,k2,k3))
   write(2, '(3F6.1,6E15.6)') kx(k1), ky(k2), kz(k3), &
   &  dble(cu1(k1,k2,k3)), imag(cu1(k1,k2,k3)), &
   &  dble(cu2(k1,k2,k3)), imag(cu2(k1,k2,k3)), &
   &  dble(cu3(k1,k2,k3)), imag(cu3(k1,k2,k3))
#endif
#ifdef field
   write(3, '(3F6.1,6E15.6)') kx(k1), ky(k2), kz(k3), &
   &  dble(a1(k1,k2,k3)), imag(a1(k1,k2,k3)), &
   &  dble(a2(k1,k2,k3)), imag(a2(k1,k2,k3)), &
   &  dble(a3(k1,k2,k3)), imag(a3(k1,k2,k3))
   write(4, '(3F6.1,6E15.6)') kx(k1), ky(k2), kz(k3), &
   &  dble(ca1(k1,k2,k3)), imag(ca1(k1,k2,k3)), &
   &  dble(ca2(k1,k2,k3)), imag(ca2(k1,k2,k3)), &
   &  dble(ca3(k1,k2,k3)), imag(ca3(k1,k2,k3))
#endif
  enddo
 enddo
enddo
#ifdef flow
close(1)
close(2)
#endif
#ifdef field
close(3)
close(4)
#endif

#ifdef flow
open(1,file='v',form='formatted')
open(2,file='cv',form='formatted')
#endif
#ifdef field
open(3,file='b',form='formatted')
open(4,file='cb',form='formatted')
#endif
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
spec=cu1
call dfftw_execute_dft_c2r(pb,spec,phys)
cv1=phys
spec=cu2
call dfftw_execute_dft_c2r(pb,spec,phys)
cv2=phys
spec=cu3
call dfftw_execute_dft_c2r(pb,spec,phys)
cv3=phys
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
spec=ca1
call dfftw_execute_dft_c2r(pb,spec,phys)
cb1=phys
spec=ca2
call dfftw_execute_dft_c2r(pb,spec,phys)
cb2=phys
spec=ca3
call dfftw_execute_dft_c2r(pb,spec,phys)
cb3=phys
#endif
! output v, b, cv and cb
do j3=1,n3
 do j2=1,n2
  do j1=1,n1
#ifdef flow
   write(1,'(6E15.6)'), x1(j1), x2(j2), x3(j3),  v1(j1,j2,j3),  v2(j1,j2,j3),  v3(j1,j2,j3)
   write(2,'(6E15.6)'), x1(j1), x2(j2), x3(j3), cv1(j1,j2,j3), cv2(j1,j2,j3), cv3(j1,j2,j3)
#endif
#ifdef field
   write(3,'(6E15.6)'), x1(j1), x2(j2), x3(j3),  b1(j1,j2,j3),  b2(j1,j2,j3),  b3(j1,j2,j3)
   write(4,'(6E15.6)'), x1(j1), x2(j2), x3(j3), cb1(j1,j2,j3), cb2(j1,j2,j3), cb3(j1,j2,j3)
#endif
  enddo
 enddo
enddo
#ifdef flow
close(1)
close(2)
#endif
#ifdef field
close(3)
close(4)
#endif
! output v in x-z plane
do j3=1,n3
 do j1=1,n1 
  write(1,'(6E15.6)') x1(j1), x3(j3), v1(j1,n2/2+1,j3), v2(j1,n2/2+1,j3), v3(j1,n2/2+1,j3), v1(j1,n2/2+1,j3)**2+v2(j1,n2/2+1,j3)**2+v3(j1,n2/2+1,j3)**2
 enddo
enddo

#ifdef field
open(1,FILE='h',form='formatted')
do j3=1,n3
 do j2=1,n2
  do j1=1,n1
   h1(j1,j2,j3)=v2(j1,j2,j3)*b3(j1,j2,j3)-v3(j1,j2,j3)*b2(j1,j2,j3)
   h2(j1,j2,j3)=v3(j1,j2,j3)*b1(j1,j2,j3)-v1(j1,j2,j3)*b3(j1,j2,j3)
   h3(j1,j2,j3)=v1(j1,j2,j3)*b2(j1,j2,j3)-v2(j1,j2,j3)*b1(j1,j2,j3)
   write(1,'(6E15.6)'), x1(j1), x2(j2), x3(j3), h1(j1,j2,j3), h2(j1,j2,j3), h3(j1,j2,j3)
  enddo
 enddo
enddo
close(1)
#endif
end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine energy(z1,z2,z3,enr1,enr2,enr3,enr)
use global
implicit none
integer j1, j2, j3, k1, k2, k3
double complex   z1(nc,n2,n3), z2(nc,n2,n3), z3(nc,n2,n3)
double precision w1(n1,n2,n3), w2(n1,n2,n3), w3(n1,n2,n3)
double precision enr1, enr2, enr3, enr
spec=z1
call dfftw_execute_dft_c2r(pb,spec,phys)
w1=phys
spec=z2
call dfftw_execute_dft_c2r(pb,spec,phys)
w2=phys
spec=z3
call dfftw_execute_dft_c2r(pb,spec,phys)
w3=phys
enr1=0.d0
enr2=0.d0
enr3=0.d0
do j3=1,n3
 do j2=1,n2
  do j1=1,n1
   enr1=enr1+w1(j1,j2,j3)**2
   enr2=enr2+w2(j1,j2,j3)**2
   enr3=enr3+w3(j1,j2,j3)**2
  enddo
 enddo
enddo
enr1=enr1/(n1*n2*n3)
enr2=enr2/(n1*n2*n3)
enr3=enr3/(n1*n2*n3)
enr=enr1+enr2+enr3
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diverg(z1,z2,z3,div)
use global
implicit none
integer j1, j2, j3, k1, k2, k3
double complex  z1(nc,n2,n3), z2(nc,n2,n3), z3(nc,n2,n3)
double precision div
div=0.d0
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   div=div+abs(kx(k1)*z1(k1,k2,k3)+ky(k2)*z2(k1,k2,k3)+kz(k3)*z3(k1,k2,k3))
  enddo
 enddo
enddo
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine curl(z1,z2,z3,cz1,cz2,cz3)
use global
implicit none
integer j1, j2, j3, k1, k2, k3
double complex  z1(nc,n2,n3),  z2(nc,n2,n3),  z3(nc,n2,n3)
double complex cz1(nc,n2,n3), cz2(nc,n2,n3), cz3(nc,n2,n3)
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   cz1(k1,k2,k3)=one*(ky(k2)*z3(k1,k2,k3)-kz(k3)*z2(k1,k2,k3))
   cz2(k1,k2,k3)=one*(kz(k3)*z1(k1,k2,k3)-kx(k1)*z3(k1,k2,k3))
   cz3(k1,k2,k3)=one*(kx(k1)*z2(k1,k2,k3)-ky(k2)*z1(k1,k2,k3))   
  enddo
 enddo
enddo
end subroutine curl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle(px1,px2,px3)
use global
implicit none
double precision px1, px2, px3
double complex vp1, vp2, vp3
integer j1, j2, j3, k1, k2, k3
vp1=0
vp2=0
vp3=0
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   vp1=vp1+u1(k1,k2,k3)*exp(one*(kx(k1)*px1+ky(k2)*px2+kz(k3)*px3))
   vp2=vp2+u2(k1,k2,k3)*exp(one*(kx(k1)*px1+ky(k2)*px2+kz(k3)*px3))
   vp3=vp3+u3(k1,k2,k3)*exp(one*(kx(k1)*px1+ky(k2)*px2+kz(k3)*px3))
  enddo
 enddo
enddo
px1=px1+dble(vp1)*dt
px2=px2+dble(vp2)*dt
px3=px3+dble(vp3)*dt
if(px1.lt.0) px1=px1+2*pi
if(px1.ge.2*pi) px1=px1-2*pi
if(px2.lt.0) px2=px2+2*pi
if(px2.ge.2*pi) px2=px2-2*pi
if(px3.lt.0) px3=px3+2*pi
if(px3.ge.2*pi) px3=px3-2*pi
end subroutine particle
