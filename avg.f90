implicit none
integer n, i, i1, i2
parameter(n=1000000)
double precision t(n),enr(n),enr1(n),enr2(n),enr3(n),ens(n),ens1(n),ens2(n),ens3(n),div
double precision enr_avg, enr1_avg, enr2_avg, enr3_avg, ens_avg, ens1_avg, ens2_avg, ens3_avg
double precision enr_dev, enr1_dev, enr2_dev, enr3_dev, ens_dev, ens1_dev, ens2_dev, ens3_dev
open(1,file='~/rotating_flow/result/5/2/flow')
do i=1,n
 read(1,'(10E15.6)') t(i),enr1(i),enr2(i),enr3(i),enr(i),ens1(i),ens2(i),ens3(i),ens(i),div
enddo
close(1)
i1=n/2
i2=n
write(6,*) 't1=',t(i1), 't2=',t(i2)
enr_avg=0
enr1_avg=0
enr2_avg=0
enr3_avg=0
ens_avg=0
ens1_avg=0
ens2_avg=0
ens3_avg=0
do i=i1,i2
 enr_avg=enr_avg+enr(i)
 enr1_avg=enr1_avg+enr1(i)
 enr2_avg=enr2_avg+enr2(i)
 enr3_avg=enr3_avg+enr3(i)
 ens_avg=ens_avg+ens(i)
 ens1_avg=ens1_avg+ens1(i)
 ens2_avg=ens2_avg+ens2(i)
 ens3_avg=ens3_avg+ens3(i)
enddo
enr_avg=enr_avg/(i2-i1+1)
enr1_avg=enr1_avg/(i2-i1+1)
enr2_avg=enr2_avg/(i2-i1+1)
enr3_avg=enr3_avg/(i2-i1+1)
write(6,*) 'averaged energy'
write(6,'(4E15.6)') enr_avg, enr1_avg, enr2_avg, enr3_avg
ens_avg=ens_avg/(i2-i1+1)
ens1_avg=ens1_avg/(i2-i1+1)
ens2_avg=ens2_avg/(i2-i1+1)
ens3_avg=ens3_avg/(i2-i1+1)
write(6,*) 'averaged dissipation'
write(6,'(4E15.6)') ens_avg, ens1_avg, ens2_avg, ens3_avg
enr_dev=0
enr1_dev=0
enr2_dev=0
enr3_dev=0
ens_dev=0
ens1_dev=0
ens2_dev=0
ens3_dev=0
do i=i1,i2
 enr_dev=enr_dev+(enr(i)-enr_avg)**2
 enr1_dev=enr1_dev+(enr1(i)-enr1_avg)**2
 enr2_dev=enr2_dev+(enr2(i)-enr2_avg)**2
 enr3_dev=enr3_dev+(enr3(i)-enr3_avg)**2
 ens_dev=ens_dev+(ens(i)-ens_avg)**2
 ens1_dev=ens1_dev+(ens1(i)-ens1_avg)**2
 ens2_dev=ens2_dev+(ens2(i)-ens2_avg)**2
 ens3_dev=ens3_dev+(ens3(i)-ens3_avg)**2
enddo
enr_dev=sqrt(enr_dev/(i2-i1+1))
enr1_dev=sqrt(enr1_dev/(i2-i1+1))
enr2_dev=sqrt(enr2_dev/(i2-i1+1))
enr3_dev=sqrt(enr3_dev/(i2-i1+1))
write(6,*) 'deviation of energy'
write(6,'(4E15.6)') enr_dev, enr1_dev, enr2_dev, enr3_dev
ens_dev=sqrt(ens_dev/(i2-i1+1))
ens1_dev=sqrt(ens1_dev/(i2-i1+1))
ens2_dev=sqrt(ens2_dev/(i2-i1+1))
ens3_dev=sqrt(ens3_dev/(i2-i1+1))
write(6,*) 'deviation of dissipation'
write(6,'(4E15.6)') ens_dev, ens1_dev, ens2_dev, ens3_dev
end
