
clear all
clc
close all
% 0.592 knot = 1 ft/s, 1 knot = 1.688 ft/s
% load factor = Total lift/ weight

d=0.0021; %density, slug/ft^3 , alt = 5000ft
cl=1.6;   %max lift coefficient, Reynolds number = 6500000
c = 4.5;  %mean aerodynamic chord, ft
b = 36.74; % wingspan
S= 175; % wing area, ft^2
AR = b^2/S; %aspect ratio 
W= 4000; %MTOW , lb
WL=W/S;  % wingload
g= 32.26; % gravity acceleration ft/s^2  
cla=4.5; % fig 3.4.14 find with AR and times 180/3.14
mu=(2*WL)/(d*g*c*cla); % mass ratio
k=(0.88*mu)/(5.3+mu); % gust alleviation factor in subsonic

u1=66; % rough air gust, ft/s, U_de
u2=50; % high speed gust, ft/s, U_de
u3=25; % dive speed gust, ft/s, U_de
vc = 193; %cruise speed, knots
vd = vc*1.25; % dive speed, typically vc*1.25
vb = vd/2; % maneuever speed at rough air gust
nlimit=2.1+24000/(W+10000); %maximum positive load factor, more than 2.1+24000/(W+10000)
nnlimit=-1.5; %maximum negative load factor

v=0 : 0.01 : vd; % velocity
n=(d*cl*((1.69*v).^2))/(2*WL); % load fator
nn=-n;

ng1 = 1+((k*u1*v*cla)/(498*WL)); %rough air gust load factor
ng2 = 1+((k*u2*v*cla)/(498*WL)); %high speed gust load factor
ng3 = 1+((k*u3*v*cla)/(498*WL)); %dive gust load factor
ng11 = 1-((k*u1*v*cla)/(498*WL)); %rough air gust load factor
ng22 = 1-((k*u2*v*cla)/(498*WL)); %high gust load factor
ng33 = 1-((k*u3*v*cla)/(498*WL)); %dive gust load factor
temp = roots([(d*cl*((1.69)^2))/(2*WL) -((k*u1*cla)/(498*WL)) -1]);
vb=temp(temp>0);

n( n > nlimit ) = nlimit; % restrict positive load factor
nn( nn < nnlimit )= nnlimit; % restrict negative load factor
n( v >= vd) = 0; % positive load factor = 0 at dive speed
nn( v >= vd) = 0; % negative load factor = 0 at dive speed

for i = find(v==vc) : find(v==vd)
    
    nn(i) = -nnlimit/(vd-vc)*(v(i)-vc)+nnlimit;
    
end

ngp1 = 1+((k*u1*vb*cla)/(498*WL)); % rough air gust point
ngp2 = 1+((k*u2*vc*cla)/(498*WL)); % high speed air gust point
ngp3 = 1+((k*u3*vd*cla)/(498*WL)); % dive speed air gust point
ngp11 = 1-((k*u1*vb*cla)/(498*WL)); % - rough air gust point 
ngp22 = 1-((k*u2*vc*cla)/(498*WL)); % - high speed air gust point
ngp33 = 1-((k*u3*vd*cla)/(498*WL)); % - dive speed air gust point

for i = 1:1001
gx(i) = 0+(i-1)*vb/1000;
gy(i) = (ngp1-1)/vb*gx(i)+1;
end
for i = 1002:2002
gx(i) = gx(1001)+(i-1002)*(vc-vb)/1000;
gy(i) = (ngp2-ngp1)/(vc-vb)*gx(i)-(ngp2-ngp1)/(vc-vb)*gx(1001)+gy(1001);
end
for i = 2003:3003
gx(i) = gx(2002)+(i-2003)*(vd-vc)/1000;
gy(i) = (ngp3-ngp2)/(vd-vc)*gx(i)-(ngp3-ngp2)/(vd-vc)*gx(2002)+gy(2002);
end

for i = 1:1001
gx2(i) = 0+(i-1)*vb/1000;
gy2(i) = (ngp11-1)/vb*gx2(i)+1;
end
for i = 1002:2002
gx2(i) = gx2(1001)+(i-1002)*(vc-vb)/1000;
gy2(i) = (ngp22-ngp11)/(vc-vb)*gx2(i)-(ngp22-ngp11)/(vc-vb)*gx2(1001)+gy2(1001);
end
for i = 2003:3003
gx2(i) = gx2(2002)+(i-2003)*(vd-vc)/1000;
gy2(i) = (ngp33-ngp22)/(vd-vc)*gx2(i)-(ngp33-ngp22)/(vd-vc)*gx2(2002)+gy2(2002);
end
gx(3004)=vd+0.001; gy(3004)=0;
gx2(3004)=vd+0.001; gy2(3004)=0;


figure(1)

% Maneuever envelope
a1 = plot(v,n, 'b','LineWidth',2);
hold on
plot(v,nn,'b', 'LineWidth',2);




% Gust envelope
 a2 = plot(gx,gy,'r', 'LineWidth',2);
 plot(gx2,gy2, 'r','LineWidth',2);
 
%Gust line 
a3 = plot(v,ng1,'m:', 'LineWidth',1);
a4 = plot(v,ng2,'m-', 'LineWidth',1);
a5 = plot(v,ng3,'m--', 'LineWidth',1);

plot(v,ng11,'m:', 'LineWidth',1)
plot(v,ng22,'m-', 'LineWidth',1)
plot(v,ng33,'m--', 'LineWidth',1)

% velocity check line
xv = vc;  
ystart = 5;
 yend = -3;
 a6 = plot([xv xv], [ystart yend],':','LineWidth',2);
 
 temp = find(n<=1);
 temp = sort(temp,'descend');
 vs = v(temp(2));
a7 = plot([vs vs], [1 -1],':','LineWidth',2);

legend([a1, a2, a3, a4, a5, a6,a7],{'Maneuever envelope', 'Gust envelope','Rough air gust',... 
        'High speed gust', 'Dive gust','Cruise speed','Stall speed'})
grid on
title('V-n diagram');
xlabel('velocity(knot)');ylabel('Load Factor');
xlim([0 vd+10])
