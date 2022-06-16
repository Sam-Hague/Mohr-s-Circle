clear
clc
clf
d1 = 60;
d2 = 30;
d3 = 45;
n1 = cos(deg2rad(d1));
n2 = cos(deg2rad(d2));
n3 = cos(deg2rad(d3));
% inputs
N_stress_x = 5;
N_stress_y = 2;
N_stress_z = 1;
S_stress_xy = 0;
S_stress_xz = 0;
S_stress_yz = 0;

% calculations
ST = [N_stress_x, S_stress_xy, S_stress_xz;
    S_stress_xy, N_stress_y, S_stress_yz;
    S_stress_xz, S_stress_yz, N_stress_z];
lambda = eig(ST);
Prin_stress_x = max(lambda);
Prin_stress_y = median(lambda);
Prin_stress_z = min(lambda);

c1 = (Prin_stress_y + Prin_stress_z)/2;
c2 = (Prin_stress_x + Prin_stress_z)/2;
c3 = (Prin_stress_x + Prin_stress_y)/2;
r1 = (abs(Prin_stress_y - Prin_stress_z))/2;
r2 = (abs(Prin_stress_x - Prin_stress_z))/2;
r3 = (abs(Prin_stress_x - Prin_stress_y))/2;
centres = [c1,c2,c3];
radii = [r1,r2,r3];
sigmatot = [Prin_stress_z, Prin_stress_y, Prin_stress_x];
% rotation point calc
sign = (Prin_stress_x*(n1*n1)) + (Prin_stress_y*(n2*n2)) + (Prin_stress_z*(n3*n3)) + ...
       2*((S_stress_xy*n1*n2)+(S_stress_xz*n1*n3)+(S_stress_yz*n2*n3));
T1 = (Prin_stress_x*n1)+(S_stress_xy*n2)+(S_stress_xz*n3);
T2 = (S_stress_xy*n1)+(Prin_stress_y*n2)+(S_stress_yz*n3);
T3 = (S_stress_xz*n1)+(S_stress_yz*n2)+(Prin_stress_z*n3);
Tn = abs(sqrt((T1*T1)+(T2*T2)+(T3*T3)-(sign*sign)));

% Min/Max
sig_min = min(sigmatot);
sig_max = max(sigmatot);
tau_max = + max(radii);
tau_min = - max(radii);

% rotation
% alpha
a = d1;
x = min(centres)-max(radii):0.01:max(centres)+max(radii);
x1 = [((tand(a)^2)+1), -2*(((tand(a)^2)*Prin_stress_z)+c1), ((Prin_stress_z^2)*(tand(a)^2)+(c1^2)-(r1^2))];
rootx1 = roots(x1);
k2 = max(rootx1);
y2 = tand(a)*(k2-Prin_stress_z);
x3 = [((tand(a)^2)+1), -2*(((tand(a)^2)*Prin_stress_z)+c2), ((Prin_stress_z^2)*(tand(a)^2)+(c2^2)-(r2^2))];
rootx3 = roots(x3);
k1 = max(rootx3);
y1 = tand(a)*(k1-Prin_stress_z);
y = tand(a)*(x-Prin_stress_z);
hold on
plot(x,y)
%Beta
b = -d2;
xb = min(centres)-max(radii):0.01:max(centres)+max(radii);
xb1 = [((tand(b)^2)+1), -2*(((tand(b)^2)*Prin_stress_x)+c3), ((Prin_stress_x^2)*(tand(b)^2)+(c3^2)-(r3^2))];
rootxb1 = roots(xb1);
kb1 = min(rootxb1);
yb1 = tand(b)*(kb1-Prin_stress_x);
xb2 = [((tand(b)^2)+1), -2*(((tand(b)^2)*Prin_stress_x)+c2), ((Prin_stress_x^2)*(tand(b)^2)+(c2^2)-(r2^2))];
rootxb2 = roots(xb2);
kb2 = min(rootxb2);
yb2 = tand(b)*(kb2-Prin_stress_x);
yb3 = tand(b)*(min(rootxb2)-Prin_stress_x);
yb = tand(b)*(xb-Prin_stress_x);
plot(xb,yb)

% charlie
c = -d3;
xc = min(centres)-max(radii):0.01:max(centres)+max(radii);
xc1 = [((tand(c)^2)+1), -2*(((tand(c)^2)*Prin_stress_y)+c1), ((Prin_stress_y^2)*(tand(c)^2)+(c1^2)-(r1^2))];
rootxc1 = roots(xc1);
kc1 = min(rootxc1);
yc1 = tand(c)*(kc1-Prin_stress_y);
xc2 = [((tand(c)^2)+1), -2*(((tand(c)^2)*Prin_stress_y)+c3), ((Prin_stress_y^2)*(tand(c)^2)+(c3^2)-(r3^2))];
rootxc2 = roots(xc2);
kc2 = max(rootxc2);
yc2 = tand(c)*(kc2-Prin_stress_y);
yc = tand(c)*(xc-Prin_stress_y);
plot(xc,yc)
% plot circles
s2 = 'g-';
th2 = 0:pi/50:2*pi;
xunit2 = r2 * cos(th2) + c2;
yunit2 = r2 * sin(th2) ;
p2 = plot(xunit2, yunit2, s2);

s1 = 'b-';
th1 = 0:pi/50:2*pi;
xunit1 = r1 * cos(th1) + c1;
yunit1 = r1 * sin(th1) ;
p1 = plot(xunit1, yunit1, s1);

s3 = 'r-';
th3 = 0:pi/50:2*pi;
xunit3 = r3 * cos(th3) + c3;
yunit3 = r3 * sin(th3) ;
p3 = plot(xunit3, yunit3, s3);

% Salient points
plot(sig_max, 0, sig_min, 0, 'Marker', 'o', 'Markersize', 6, 'color','Black')
plot(c2, tau_max, c2, tau_min, 'Marker', 'o', 'Markersize', 6, 'color', 'Black')
 plot(k1,y1,k2,y2,c3,0,'Marker', 'o', 'Markersize', 6, 'color', 'Red')
 plot(kb2,yb2,kb1,yb1,c1,0,'Marker','+','Markersize', 6,'color','blue')
 plot(kc2,yc2,kc1,yc1,c2,0,'Marker','*','Markersize', 6,'color','yellow')
plot(sign,Tn,'Marker','o','Markersize', 10, 'color', 'Black')
axis([min(centres)-max(radii),max(centres)+max(radii), tau_min-(max(radii))/10, tau_max+(max(radii))/10])

% arc drawing alpha
xa = 100.3343;
ya = 0;
ra = sqrt(((100.3343-1.6)^2)+((0-38.57)^2));
tha = 0:pi/1000:2*pi;
xunit = ra * cos(tha) + xa;
xshort = (xunit(xunit>=k2 & xunit<=k1));
xhalf = xshort(1,1:(round(length(xshort))/2));
yunit = ra * sin(tha) + ya;
yshort = (yunit(yunit>=y2 & yunit<=y1));
yhalf = yshort(1,(round(length(yshort)/2))+1:end);

% arc drawing Beta
xb = c1;
yb = 0;
rb = sqrt(((xb-kb1)^2+(0-yb1)^2));
thb = 0:pi/1000:2*pi;
xunitb = rb * cos(thb) + xb;
xshortb = (xunitb(xunitb>=kb2 & xunitb<kb1));
xhalfb = xshortb(1,1:(round(length(xshortb))/2));
yunitb = rb * sin(thb) + yb;
yshortb = (yunitb(yunitb>=yb1 & yunitb<=yb2));
yhalfb = flip(yshortb(1,(round(length(yshortb)/2))+1:end));

% arc drawing Charlie
th1 = tan(yc1/(c2+kc1)); %0.625628...
th2 = pi+tan(yc2/(kc2-c2)); %2.2699...
xc = c2;
yc = 0;
rc = sqrt(((xc-kc1)^2+(0-yc1)^2));
thc = -th1:pi/100:th2;
xunitc = rc * cos(thc) + xc;
yunitc = rc * sin(thc) + yc;
plot(xunitc,yunitc)

% plotting
plot(xhalf, yhalf, 'color', 'black')
plot(xhalfb, yhalfb, 'color', 'blue')

plot([k1,xa,k2], [y1,ya,y2], 'black:')
plot([kb1,xb,kb2], [yb1,yb,yb2],'blue:')
plot([kc1,xc,kc2], [yc1,yc,yc2],'green:')
%plot(xshort,yshort)
%plot(xunit, yunit)