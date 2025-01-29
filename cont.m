clear; close all; clc;

n1=64; n3=64;

load xz
xz1=xz;
x1=xz1(:,1);
y1=xz1(:,2);
z1=xz1(:,3);
x1=reshape(x1,n1,n3);
y1=reshape(y1,n1,n3);
z1=reshape(z1,n1,n3);

%load xz4
%x4=xz4(:,1);
%y4=xz4(:,2);
%z4=xz4(:,6);
%x4=reshape(x4,n1,n3);
%y4=reshape(y4,n1,n3);
%z4=reshape(z4,n1,n3);

x0=[0 2*pi 2*pi 0 0];
y0=[0 0 2*pi 2*pi 0];

figure
axis equal
axis off
hold

xc=0; yc=0;
plot(x0+xc,y0+yc,'-k')
contour(x1+xc,y1+yc,z1,'-k')

%xc=2*pi+1; yc=0;
%plot(x0+xc,y0+yc,'-k')
%contour(x4+xc,y4+yc,z4,'-k')
