! wavenumber 0, 1 ... n/2-1, n/2, -(n/2-1) ... -1 corresponds to index 1, 2 ... n/2, n/2+1, n/2+2 ... n of enr and ens

#define flow
#define field

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

! velocity variables
!double complex    u1(nc,n2,n3),  u2(nc,n2,n3),  u3(nc,n2,n3) ! velocity in spectral space
!double precision  v1(n1,n2,n3),  v2(n1,n2,n3),  v3(n1,n2,n3) ! velocity in physical space
!double complex   cu1(nc,n2,n3), cu2(nc,n2,n3), cu3(nc,n2,n3) ! vorticity in spectral space
!double precision cv1(n1,n2,n3), cv2(n1,n2,n3), cv3(n1,n2,n3) ! vorticity in physical space

! field variables
!double complex    a1(nc,n2,n3),  a2(nc,n2,n3),  a3(nc,n2,n3) ! field in spectral space
!double precision  b1(n1,n2,n3),  b2(n1,n2,n3),  b3(n1,n2,n3) ! field in physical space
!double complex   ca1(nc,n2,n3), ca2(nc,n2,n3), ca3(nc,n2,n3) ! current in spectral space
!double precision cb1(n1,n2,n3), cb2(n1,n2,n3), cb3(n1,n2,n3) ! current in physical space

integer j1, j2, j3, k1, k2, k3
double precision i1, i2, i3, w1r, w1i, w2r, w2i, w3r, w3i
double precision enr(n1,n2,n3), ens(n1,n2,n3)
double precision enrsum, enrmax, enssum, ensmax
integer kx, ky, kz
double precision x, y

double precision enr1(n1/2+1,n2/2+1,n3/2+1)
double precision ens1(n1/2+1,n2/2+1,n3/2+1)
integer wavenum2(0:n1/2,0:n2/2,0:n3/2)
double precision enr2(0:(n1/2)**2+(n2/2)**2+(n3/2)**2)
double precision ens2(0:(n1/2)**2+(n2/2)**2+(n3/2)**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef flow

! read in flow
open(1,file='u',form='formatted')
open(2,file='cu',form='formatted')
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   read(1, '(3F6.1,6E15.6)') i1, i2, i3, w1r, w1i, w2r, w2i, w3r, w3i
   enr(k1,k2,k3)=(w1r*w1r+w1i*w1i+w2r*w2r+w2i*w2i+w3r*w3r+w3i*w3i)
   read(2, '(3F6.1,6E15.6)') i1, i2, i3, w1r, w1i, w2r, w2i, w3r, w3i
   ens(k1,k2,k3)=(w1r*w1r+w1i*w1i+w2r*w2r+w2i*w2i+w3r*w3r+w3i*w3i)
  enddo
 enddo
enddo
close(1)
close(2)

! give the other half
do k3=1,n3
 do k2=1,n2
  do k1=nc+1,n1
   enr(k1,k2,k3)=enr(n1-k1+2,k2,k3)
   ens(k1,k2,k3)=ens(n1-k1+2,k2,k3)
  enddo
 enddo
enddo

! summation and maximum
enrsum=0.d0
enrmax=0.d0
do k3=1,n3
 do k2=1,n2
  do k1=1,n1
   enrsum=enrsum+enr(k1,k2,k3)
   if(enr(k1,k2,k3).ge.enrmax) then
    enrmax=enr(k1,k2,k3)
    kx=k1-1
    if((k1-1).gt.n1/2) kx=k1-1-n1
    ky=k2-1
    if((k2-1).gt.n2/2) ky=k2-1-n2
    kz=k3-1
    if((k3-1).gt.n3/2) kz=k3-1-n3
   endif
  enddo
 enddo
enddo
write(6,*) 'sum of u^2=',enrsum,'max=',enrmax,'on','(',kx, ky, kz,')'

enssum=0.d0
ensmax=0.d0
do k3=1,n3
 do k2=1,n2
  do k1=1,n1
   enssum=enssum+ens(k1,k2,k3)
   if(ens(k1,k2,k3).ge.ensmax) then
    ensmax=ens(k1,k2,k3)
    kx=k1-1
    if((k1-1).gt.n1/2) kx=k1-1-n1
    ky=k2-1
    if((k2-1).gt.n2/2) ky=k2-1-n2
    kz=k3-1
    if((k3-1).gt.n3/2) kz=k3-1-n3
   endif
  enddo
 enddo
enddo
write(6,*), 'sum of w^2=',enssum,'max=',ensmax,'on','(',kx, ky, kz,')'

! u^2 on kx
enrsum=0
do k1=1,n1/2+1
 x=k1-1
 y=0
 do k2=1,n2
  do k3=1,n3
   y=y+enr(k1,k2,k3)
   if((k1.ne.1).and.(k1.ne.n1/2+1)) y=y+enr(n1-k1+2,k2,k3)
  enddo
 enddo
 write(11,*) x, y
 enrsum=enrsum+y
enddo
write(6,*) 'sum of u^2 on kx', enrsum

! w^2 on kx
enssum=0
do k1=1,n1/2+1
 x=k1-1
 y=0
 do k2=1,n2
  do k3=1,n3
   y=y+ens(k1,k2,k3)
   if((k1.ne.1).and.(k1.ne.n1/2+1)) y=y+ens(n1-k1+2,k2,k3)
  enddo
 enddo
 write(14,*) x, y
 enssum=enssum+y
enddo
write(6,*) 'sum of w^2 on kx', enssum

! u^2 on ky
enrsum=0
do k2=1,n2/2+1
 x=k2-1
 y=0
 do k1=1,n1
  do k3=1,n3
   y=y+enr(k1,k2,k3)
   if((k2.ne.1).and.(k2.ne.n2/2+1)) y=y+enr(k1,n2-k2+2,k3)
  enddo
 enddo
 write(12,*) x, y
 enrsum=enrsum+y
enddo
write(6,*) 'sum of u^2 on ky', enrsum

! w^2 on ky
enssum=0
do k2=1,n2/2+1
 x=k2-1
 y=0
 do k1=1,n1
  do k3=1,n3
   y=y+ens(k1,k2,k3)
   if((k2.ne.1).and.(k2.ne.n2/2+1)) y=y+ens(k1,n2-k2+2,k3)
  enddo
 enddo
 write(15,*) x, y
 enssum=enssum+y
enddo
write(6,*) 'sum of w^2 on ky', enssum

! u^2 on kz
enrsum=0
do k3=1,n3/2+1
 x=k3-1
 y=0
 do k1=1,n1
  do k2=1,n2
   y=y+enr(k1,k2,k3)
   if((k3.ne.1).and.(k3.ne.n3/2+1)) y=y+enr(k1,k2,n3-k3+2)
  enddo
 enddo
 write(13,*) x, y
 enrsum=enrsum+y
enddo
write(6,*) 'sum of u^2 on kz', enrsum

! w^2 on kz
enssum=0
do k3=1,n3/2+1
 x=k3-1
 y=0
 do k1=1,n1
  do k2=1,n2
   y=y+ens(k1,k2,k3)
   if((k3.ne.1).and.(k3.ne.n3/2+1)) y=y+ens(k1,k2,n3-k3+2)
  enddo
 enddo
 write(16,*) x, y
 enssum=enssum+y
enddo
write(6,*) 'sum of w^2 on kz', enssum

! in wavenumber space
do k1=1,n1/2+1
 do k2=1,n2/2+1
  do k3=1,n3/2+1
   enr1(k1,k2,k3)=enr(k1,k2,k3)
   ens1(k1,k2,k3)=ens(k1,k2,k3)
   if((k1.ne.1).and.(k1.ne.n1/2+1)) enr1(k1,k2,k3)=enr1(k1,k2,k3)+enr(n1-k1+2,k2,k3)
   if((k2.ne.1).and.(k2.ne.n2/2+1)) enr1(k1,k2,k3)=enr1(k1,k2,k3)+enr(k1,n2-k2+2,k3)
   if((k3.ne.1).and.(k3.ne.n3/2+1)) enr1(k1,k2,k3)=enr1(k1,k2,k3)+enr(k1,k2,n3-k3+2)
   if((k1.ne.1).and.(k1.ne.n1/2+1)) ens1(k1,k2,k3)=ens1(k1,k2,k3)+ens(n1-k1+2,k2,k3)
   if((k2.ne.1).and.(k2.ne.n2/2+1)) ens1(k1,k2,k3)=ens1(k1,k2,k3)+ens(k1,n2-k2+2,k3)
   if((k3.ne.1).and.(k3.ne.n3/2+1)) ens1(k1,k2,k3)=ens1(k1,k2,k3)+ens(k1,k2,n3-k3+2)
  enddo
 enddo
enddo

do k1=0,n1/2
 do k2=0,n2/2
  do k3=0,n3/2
   wavenum2(k1,k2,k3)=k1**2+k2**2+k3**2
  enddo
 enddo
enddo

do j1=0,(n1/2)**2+(n2/2)**2+(n3/2)**2
 enr2(j1)=0.d0
 ens2(j1)=0.d0
 do k1=0,n1/2
  do k2=0,n2/2
   do k3=0,n3/2
    if(wavenum2(k1,k2,k3).eq.j1) then
     enr2(j1)=enr2(j1)+enr1(k1+1,k2+1,k3+1)
     ens2(j1)=ens2(j1)+ens1(k1+1,k2+1,k3+1)
    endif
   enddo
  enddo
 enddo
 write(10,'(3E25.15)') sqrt(dble(j1)), enr2(j1), ens2(j1)
enddo
close(10)
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef field

! read in field
open(3,file='a',form='formatted')
open(4,file='ca',form='formatted')
do k3=1,n3
 do k2=1,n2
  do k1=1,nc
   read(3, '(3F6.1,6E15.6)') i1, i2, i3, w1r, w1i, w2r, w2i, w3r, w3i
   enr(k1,k2,k3)=(w1r*w1r+w1i*w1i+w2r*w2r+w2i*w2i+w3r*w3r+w3i*w3i)/2
   read(4, '(3F6.1,6E15.6)') i1, i2, i3, w1r, w1i, w2r, w2i, w3r, w3i
   ens(k1,k2,k3)=(w1r*w1r+w1i*w1i+w2r*w2r+w2i*w2i+w3r*w3r+w3i*w3i)/2
  enddo
 enddo
enddo
close(3)
close(4)

! give the other half
do k3=1,n3
 do k2=1,n2
  do k1=nc+1,n1
   enr(k1,k2,k3)=enr(n1-k1+2,k2,k3)
   ens(k1,k2,k3)=ens(n1-k1+2,k2,k3)
  enddo
 enddo
enddo

! summation and maximum
enrsum=0
enrmax=0
do k3=1,n3
 do k2=1,n2
  do k1=1,n1
   enrsum=enrsum+enr(k1,k2,k3)
   if(enr(k1,k2,k3).ge.enrmax) then
    enrmax=enr(k1,k2,k3)
    kx=k1-1
    if((k1-1).gt.n1/2) kx=k1-1-n1
    ky=k2-1
    if((k2-1).gt.n2/2) ky=k2-1-n2
    kz=k3-1
    if((k3-1).gt.n3/2) kz=k3-1-n3
   endif
  enddo
 enddo
enddo
write(6,*), 'sum of b^2=',enrsum,'max=',enrmax,'on','(',kx, ky, kz,')'

enssum=0.d0
ensmax=0.d0
do k3=1,n3
 do k2=1,n2
  do k1=1,n1
   enssum=enssum+ens(k1,k2,k3)
   if(ens(k1,k2,k3).ge.ensmax) then
    ensmax=ens(k1,k2,k3)
    kx=k1-1
    if((k1-1).gt.n1/2) kx=k1-1-n1
    ky=k2-1
    if((k2-1).gt.n2/2) ky=k2-1-n2
    kz=k3-1
    if((k3-1).gt.n3/2) kz=k3-1-n3
   endif
  enddo
 enddo
enddo
write(6,*), 'sum of j^2=',enssum,'max=',ensmax,'on','(',kx, ky, kz,')'

! b^2 on kx
enrsum=0
do k1=1,n1/2+1
 x=k1-1
 y=0
 do k2=1,n2
  do k3=1,n3
   y=y+enr(k1,k2,k3)
   if((k1.ne.1).and.(k1.ne.n1/2+1)) y=y+enr(n1-k1+2,k2,k3)
  enddo
 enddo
 write(21,*) x, y
 enrsum=enrsum+y
enddo
write(6,*) 'sum of b^2 on kx', enrsum

! j^2 on kx
enssum=0
do k1=1,n1/2+1
 x=k1-1
 y=0
 do k2=1,n2
  do k3=1,n3
   y=y+ens(k1,k2,k3)
   if((k1.ne.1).and.(k1.ne.n1/2+1)) y=y+ens(n1-k1+2,k2,k3)
  enddo
 enddo
 write(24,*) x, y
 enssum=enssum+y
enddo
write(6,*) 'sum of j^2 on kx', enssum

! b^2 on ky
enrsum=0
do k2=1,n2/2+1
 x=k2-1
 y=0
 do k1=1,n1
  do k3=1,n3
   y=y+enr(k1,k2,k3)
   if((k2.ne.1).and.(k2.ne.n2/2+1)) y=y+enr(k1,n2-k2+2,k3)
  enddo
 enddo
 write(22,*) x, y
 enrsum=enrsum+y
enddo
write(6,*) 'sum of b^2 on ky', enrsum

! j^2 on ky
enssum=0
do k2=1,n2/2+1
 x=k2-1
 y=0
 do k1=1,n1
  do k3=1,n3
   y=y+ens(k1,k2,k3)
   if((k2.ne.1).and.(k2.ne.n2/2+1)) y=y+ens(k1,n2-k2+2,k3)
  enddo
 enddo
 write(25,*) x, y
 enssum=enssum+y
enddo
write(6,*) 'sum of j^2 on ky', enssum

! b^2 on kz
enrsum=0
do k3=1,n3/2+1
 x=k3-1
 y=0
 do k1=1,n1
  do k2=1,n2
   y=y+enr(k1,k2,k3)
   if((k3.ne.1).and.(k3.ne.n3/2+1)) y=y+enr(k1,k2,n3-k3+2)
  enddo
 enddo
 write(23,*) x, y
 enrsum=enrsum+y
enddo
write(6,*) 'sum of b^2 on kz', enrsum

! j^2 on kz
enssum=0
do k3=1,n3/2+1
 x=k3-1
 y=0
 do k1=1,n1
  do k2=1,n2
   y=y+ens(k1,k2,k3)
   if((k3.ne.1).and.(k3.ne.n3/2+1)) y=y+ens(k1,k2,n3-k3+2)
  enddo
 enddo
 write(26,*) x, y
 enssum=enssum+y
enddo
write(6,*) 'sum of j^2 on kz', enssum

! in wavenumber space
do k1=1,n1/2+1
 do k2=1,n2/2+1
  do k3=1,n3/2+1
   enr1(k1,k2,k3)=enr(k1,k2,k3)
   ens1(k1,k2,k3)=ens(k1,k2,k3)
   if((k1.ne.1).and.(k1.ne.n1/2+1)) enr1(k1,k2,k3)=enr1(k1,k2,k3)+enr(n1-k1+2,k2,k3)
   if((k2.ne.1).and.(k2.ne.n2/2+1)) enr1(k1,k2,k3)=enr1(k1,k2,k3)+enr(k1,n2-k2+2,k3)
   if((k3.ne.1).and.(k3.ne.n3/2+1)) enr1(k1,k2,k3)=enr1(k1,k2,k3)+enr(k1,k2,n3-k3+2)
   if((k1.ne.1).and.(k1.ne.n1/2+1)) ens1(k1,k2,k3)=ens1(k1,k2,k3)+ens(n1-k1+2,k2,k3)
   if((k2.ne.1).and.(k2.ne.n2/2+1)) ens1(k1,k2,k3)=ens1(k1,k2,k3)+ens(k1,n2-k2+2,k3)
   if((k3.ne.1).and.(k3.ne.n3/2+1)) ens1(k1,k2,k3)=ens1(k1,k2,k3)+ens(k1,k2,n3-k3+2)
  enddo
 enddo
enddo

do k1=0,n1/2
 do k2=0,n2/2
  do k3=0,n3/2
   wavenum2(k1,k2,k3)=k1**2+k2**2+k3**2
  enddo
 enddo
enddo

do j1=0,(n1/2)**2+(n2/2)**2+(n3/2)**2
 enr2(j1)=0.d0
 ens2(j1)=0.d0
 do k1=0,n1/2
  do k2=0,n2/2
   do k3=0,n3/2
    if(wavenum2(k1,k2,k3).eq.j1) then
     enr2(j1)=enr2(j1)+enr1(k1+1,k2+1,k3+1)
     ens2(j1)=ens2(j1)+ens1(k1+1,k2+1,k3+1)
    endif
   enddo
  enddo
 enddo
 write(20,'(3E25.15)') sqrt(dble(j1)), enr2(j1), ens2(j1)
enddo
close(20)
close(21)
close(22)
close(23)
close(24)
close(25)
close(26)
#endif

end program main
