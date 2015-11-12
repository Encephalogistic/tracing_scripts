function [sl,az,error,stdr,B,sllb,slub,azlb,azub, BINT, R, RINT, STATS]=regressplanecone(x,y,z)
%Modified from regressplane.m, by Kevin
%regress solves yt. Here, 'y' is z and 'X' is made of the three columns
%of the x vector, the y vector, and a column of ones for the constant term
X=[x y ones(length(y),1)];
[B,BINT,R,RINT,STATS] = regress(z,X,.05);

sl=atand((B(1)^2+B(2)^2)/(sqrt(B(1)^2+B(2)^2)));
az=atan2d(-B(2),-B(1));

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

% Cone Begins

slub = sl;
sllb = sl;
azub = az;
azlb = az;

n = cross([1,0,B(1)],[0,1,B(2)]);
r = norm(n)*sind(error);
nu = n/norm(n);
if n(1) ~= 1 && n(1) ~= -1
    u = cross(n,[1,0,0])/norm(cross(n,[1,0,0]));
else
    u = cross(n,[0,1,0])/norm(cross(n,[0,1,0]));
end

[a,b] = meshgrid([-max(n(1))*5 max(n(1))*5],[-max(n(2))*5 max(n(2))*5]);

%figure
%hold on
%surf(a,b,B(1)*a + B(2)*b + B(3),'EdgeColor','none')
%quiver3(0,0,B(3),n(1),n(2),n(3))
%quiver3(n(1),n(2),n(3)+B(3),u(1),u(2),u(3))

%figure
%hold on

for i = 1:360
    p = r*cosd(i)*u + r*sind(i)*cross(nu,u) + n;
    slub = max(slub,atand((p(1)^2+p(2)^2)/sqrt(p(1)^2+p(2)^2)));
    sllb = min(sllb,atand((p(1)^2+p(2)^2)/sqrt(p(1)^2+p(2)^2)));
    if az > 0
        azub = min(azub,atan2d(-p(2),-p(1))+180);
        azlb = min(azlb,atan2d(-p(2),-p(1))+180);
        temp = atan2d(-p(2),-p(1))+180;
    else
        azub = min(azub,atan2d(-p(2),-p(1))-180);
        azlb = min(azlb,atan2d(-p(2),-p(1))-180);
        temp = atan2d(-p(2),-p(1))-180;
    end
    %scatter(atand((p(1)^2+p(2)^2)/sqrt(p(1)^2+p(2)^2)),i)
    %scatter3(p(1),p(2),p(3)+B(3))
end


