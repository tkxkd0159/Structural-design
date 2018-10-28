
clear all
clc
close all
% 1 knot = 1.688 ft/s
d=0.00089068; %density, slug/ft^3 , alt = 30000ft
cl=1.3;   %lift coefficient
c = 8.53;  %mean aerodynamic chord
S= 890; % wing area, ft^2
W= 66000; %MTOW , lb
WL=66000/890;  % wingload
g= 32.26; % gravity acceleration ft/s^2  
a=2; % CL_alpha
mu=(2*WL)/(d*g*c*a); % mass ratio
k=(0.88*mu)/(5.3+mu); % gust alleviation factor in subsonic

u1=66; % rough air gust, ft/s, U_de
u2=50; % high speed gust, ft/s, U_de
u3=25; % dive speed gust, ft/s, U_de
vc=640; %cruise speed, knots
vd=vc*1.25; % dive speed
nlimit=2.5;
ndivelimit=-1;
velocity=0 : 0.1 : vd + 1.1;
n=(d*cl*(velocity.^2))/(2*WL);
n2=(d*cl*(velocity.^2))/(2*WL);
ndive=-n;

ng1 = 1+((k*u1*velocity*a)/(498*WL));
ng2 = 1+((k*u2*velocity*a)/(498*WL));
ng3 = 1+((k*u3*velocity*a)/(498*WL));
ng11 = 1-((k*u1*velocity*a)/(498*WL));
ng22 = 1-((k*u2*velocity*a)/(498*WL));
ng33 = 1-((k*u3*velocity*a)/(498*WL));

inx = (velocity <= vd);
n( n>= nlimit ) = nlimit;
n( velocity >= vd) = 0;
ndive( velocity >= vd) = 0;
ndive( ndive <= ndivelimit )= ndivelimit;
tmp = ndivelimit+((-ndivelimit./(vd-vc)).*(velocity-vc));
ndive( velocity >= vc) = tmp(velocity >= vc);
na=(d*cl*(velocity.^2))/(2*WL);
ndivea=-na;
inxa = (velocity <= vd);
na( na>=2.828) = 2.828;
tm = 3.169+(((3.169-2.828)./(380-242.6)).*(velocity-380));
na( velocity >=242.6) = tm(velocity >= 242.6);
mp = 2.356+(((2.356-3.169)./(vd-vc)).*(velocity-vd));
na( velocity >=380) = mp(velocity >= 380);
na( velocity >= vd) = 0;
ndivea( ndivea <= ndivelimit )= ndivelimit;
tmp1 = (-1)+(((1-1.169)./(vc-265.5)).*(velocity-265.5));
ndivea( velocity >= 265.5) = tmp1(velocity >= 265.5);
tmp = -1.169+(((1.169-0.3556)./(vd-vc)).*(velocity-vc));
ndivea( velocity >= vc) = tmp(velocity >= vc);
ndivea( velocity >= vd) = 0;

figure(1);
plot(velocity(inxa),ndivea(inxa),'r', 'LineWidth',2)
hold on
plot(velocity(inxa),na(inxa), 'r','LineWidth',2);
plot(velocity(inx),ndive(inx),'b', 'LineWidth',2)
hold on
plot(velocity(inx),n(inx), 'b','LineWidth',2);
hold on;
plot(velocity,ng1,'r:', 'LineWidth',2)
hold on;
plot(velocity,ng2,'r:', 'LineWidth',2)
hold on;
plot(velocity,ng3,'r:', 'LineWidth',2)
hold on;
grid on;
title('V-n diagram');
xlabel('velocity');ylabel('Load Factor');
