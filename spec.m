close all
clear
clc

n1=64;
n2=64;
n3=64;

load fort.11
u2kx=fort;
load fort.12
u2ky=fort;
load fort.13
u2kz=fort;
load fort.14
w2kx=fort;
load fort.15
w2ky=fort;
load fort.16
w2kz=fort;

figure
xmin=0;x1max=n1/2;x2max=n2/2;x3max=n3/2;
subplot(2,3,1)
plot(u2kx(:,1),u2kx(:,2),'-k')
xlim([xmin x1max]);
subplot(2,3,2)
plot(u2ky(:,1),u2ky(:,2),'-k')
xlim([xmin x2max]);
subplot(2,3,3)
plot(u2kz(:,1),u2kz(:,2),'-k')
xlim([xmin x3max]);
subplot(2,3,4)
plot(w2kx(:,1),w2kx(:,2),'-k')
xlim([xmin x1max]);
subplot(2,3,5)
plot(w2ky(:,1),w2ky(:,2),'-k')
xlim([xmin x2max]);
subplot(2,3,6)
plot(w2kz(:,1),w2kz(:,2),'-k')
xlim([xmin x3max]);

load fort.21
b2kx=fort;
load fort.22
b2ky=fort;
load fort.23
b2kz=fort;
load fort.24
j2kx=fort;
load fort.25
j2ky=fort;
load fort.26
j2kz=fort;

figure
xmin=0;x1max=n1/2;x2max=n2/2;x3max=n3/2;
subplot(2,3,1)
plot(b2kx(:,1),b2kx(:,2),'-k')
xlim([xmin x1max]);
subplot(2,3,2)
plot(b2ky(:,1),b2ky(:,2),'-k')
xlim([xmin x2max]);
subplot(2,3,3)
plot(b2kz(:,1),b2kz(:,2),'-k')
xlim([xmin x3max]);
subplot(2,3,4)
plot(j2kx(:,1),j2kx(:,2),'-k')
xlim([xmin x1max]);
subplot(2,3,5)
plot(j2ky(:,1),j2ky(:,2),'-k')
xlim([xmin x2max]);
subplot(2,3,6)
plot(j2kz(:,1),j2kz(:,2),'-k')
xlim([xmin x3max]);
