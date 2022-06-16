clear
clc
Theta = 30;
cx1 = 2;
cy1 = 0;
cx2 = 0;
cy2 = 0;
cx3 = 6;
cy3 = 0;
r1 = 6;
r2 = 4;
r3 = 2;
th = 0: pi/1000 : 2*pi;
x1 = r1*sin(th)+cx1;
y1 = r1*cos(th)+cy1;
x2 = r2*sin(th)+cx2;
y2 = r2*cos(th)+cy2;
x3 = r3*sin(th)+cx3;
y3 = r3*cos(th)+cy3;
N_stress_x = 30;
N_stress_y = 40;
S_stress_yface = 10;
cx = (N_stress_y+N_stress_x)/2;
cy = 0;
r = sqrt(((N_stress_x-cx)^2)+(S_stress_yface^2));
ThetaI = atand((S_stress_yface)/(N_stress_x-cx));
% rotated element
ThetaT = (Theta)+ThetaI;
ya = -r*sind(ThetaT)+cy;
xa = r*cosd(ThetaT)+cx;
yb = -ya;
xb = cx-(xa-cx);
arcr = 10;
arcth = deg2rad(ThetaI):pi/100:deg2rad(ThetaT);
[xx, yy] = pol2cart(arcth,arcr);
xx(2:end+1) = xx(1:end);
xx(1) = 0;
xx(end) = 0;
yy(2:end+1) = yy(1:end);
yy(1) = 0;
yy(end) = 0;
[tt, hh] = pol2cart(arcth(length(arcth/2)),arcr+5);


arc1x = r1*sin(0:pi/100:2*pi)+cx1;
arc1y = r1*cos(0:pi/100:2*pi)+cy1;
arc2x = r2*sin(0:pi/100:2*pi)+cx2;
arc2y = r2*cos(0:pi/100:2*pi)+cy2;
arc3x = r3*sin(0:pi/100:2*pi)+cx3;
arc3y = r3*cos(0:pi/100:2*pi);
hold on
plot(arc1x,arc1y)
plot(arc2x,arc2y)
plot(arc3x,arc3y)
patch(xx+cx, -yy+cy,'g','edgecolor','none')
xtot = [arc2x arc3x];
ytot = [arc2y arc3y];
fill(arc1x, arc1y, 'r')
fill(xtot,ytot,'w')
axis("equal")
