function [sl,az,error,stdr,B,sllb,slub,azlb,azub, BINT, R, RINT, STATS]=regressplanemod(x,y,z)
%Modified from regressplane.m, by Kevin
%regress solves yt. Here, 'y' is z and 'X' is made of the three columns
%of the x vector, the y vector, and a column of ones for the constant term
X=[x y ones(length(y),1)];
[B,BINT,R,RINT,STATS] = regress(z,X,.05);

sl=atand((B(1)^2+B(2)^2)/(sqrt(B(1)^2+B(2)^2)));
slub=atand((BINT(1,1)^2+BINT(2,1)^2)/(sqrt(BINT(1,1)^2+BINT(2,1)^2)));
sllb=atand((BINT(1,2)^2+BINT(2,2)^2)/(sqrt(BINT(1,2)^2+BINT(2,2)^2)));
az=atan2d(-B(2),-B(1));
azlb=atan2d(-BINT(2,1),-BINT(1,1));
azub=atan2d(-BINT(2,2),-BINT(1,2));

error(1)=acosd(dot([B(1),B(2),-1],[BINT(1,1) BINT(2,1) -1])/(norm([B(1),B(2), ...
		    -1])*norm([BINT(1,1) BINT(2,1) -1])));
error(2)=acosd(dot([B(1),B(2),-1],[BINT(1,1) BINT(2,2) -1])/(norm([B(1),B(2), ...
		    -1])*norm([BINT(1,1) BINT(2,2) -1])));
error(3)=acosd(dot([B(1),B(2),-1],[BINT(1,2) BINT(2,1) -1])/(norm([B(1),B(2), ...
		    -1])*norm([BINT(1,2) BINT(2,1) -1])));
error(4)=acosd(dot([B(1),B(2),-1],[BINT(1,2) BINT(2,2) -1])/(norm([B(1),B(2), ...
		    -1])*norm([BINT(1,2) BINT(2,2) -1])));
error=max(error);

stdr=std(R);