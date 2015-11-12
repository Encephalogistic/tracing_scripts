function [sl,az,error,stdr,B]=regressplane(x,y,z)
%regress solves yt. Here, 'y' is z and 'X' is made of the three columns
%of the x vector, the y vector, and a column of ones for the constant term
X=[x y ones(length(y),1)];
[B,BINT,R,RINT,STATS] = regress(z,X,.05);
%sl=180/pi*atan((abs(B(1))+abs(B(2)))/sqrt(2));
%sl=180/pi*atan(sqrt(B(1).^2+B(2).^2));


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
%old code, more correct to just calculate error of normal vector,
%particularly for instances where error region includes horizontal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%errsl1=abs(sl-180/pi*atan((BINT(1,1)^2+BINT(2,1)^2)/(sqrt(BINT(1,1)^2+BINT(2,1)^2))));
%errsl2=abs(sl-180/pi*atan((BINT(1,2)^2+BINT(2,2)^2)/(sqrt(BINT(1,2)^2+BINT(2,2)^2))));

%dsl=max([errsl1,errsl2]);

%erraz1=abs(180/pi*atan2(-BINT(2,1),-BINT(1,1))-az);
%if erraz1>180, erraz1=360-erraz1; end

%erraz2=abs(180/pi*atan2(-BINT(2,2),-BINT(1,2))-az);
%if erraz2>180, erraz2=360-erraz2; end

%daz=max([erraz1,erraz2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dsl=1/(1+B(1)^2+B(2)^2)*(2^-.5)*(abs((BINT(1,1)-B(1))/B(1))+abs((BINT(2,1)-B(2))/B(2)))*(B(1)^2+B(2)^2)^.5*180/pi;
%daz=1/(1+(B(2)/B(1))^2)*(abs((BINT(1,1)-B(1))/B(1))^2+abs((BINT(2,1)-B(2))/B(2))^2)^.5*(B(2)/B(1))*180/pi;
%km to meters conversion...
stdr=std(R);
