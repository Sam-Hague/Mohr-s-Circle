
xa = 100.3343;
ya = 0;
ra = sqrt(((100.3343-1.6)^2)+((0-38.57)^2));

tha = 0:pi/50:pi;
xunit = ra * cos(tha) + xa;
yunit = ra * sin(tha) + ya;
plot(xunit,yunit)