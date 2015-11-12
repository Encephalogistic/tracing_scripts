%function [relev, angle, dist, slope, azimuth, xmeters, ymeters] = traceToPolar(fileName,container,mound,mat,r,a)

%============
% File inputs:
%   fileName    => the name of the shapefile group containing layer trace
%   container   => 1-40, the container surrounding the trace
%   mound       => 1-n, identifies which mound in the container to use as
%                   origin
%============

fileName = 'GCLT_1750_VTP';
container = 8;
mound = 1;

%============

mx = [];
my = [];
ras = [];
maxId = 0;
divS = struct('mx',mx,'my',my,'ras',ras);

newShape = shaperead(fileName);

for i = 1:length(newShape)
    maxId = max(maxId,newShape(i).Id);
end

stretchx = cos(ccent(container).cc(2)/3376200.00386);

divS(maxId+1).mx = [];
slope(maxId+1) = 0;
azimuth(maxId+1) = 0;
error(maxId+1) = 0;
stdr(maxId+1) = 0;
BVec(maxId+1).B = 0;

xmeters(maxId+1) = 0;
ymeters(maxId+1) = 0;
relev(maxId+1) = 0;
angle(maxId+1) = 0;
dist(maxId+1) = 0;
distc(maxId+1) = 0;

for i=1:maxId+1
    distc(i) = nan;
end

%read in old data
if exist(['polarData_C',num2str(container),'M',num2str(mound),'.csv'],'file')
    allFiles = importdata(['polarData_C',num2str(container),'M',num2str(mound),'.csv'],',');
else
    csvwrite(['polarData_C',num2str(container),'M',num2str(mound),'.csv'],[]);
    ['File created: polarData_C',num2str(container),'M',num2str(mound),'.csv']
    allFiles = importdata(['polarData_C',num2str(container),'M',num2str(mound),mound,'.csv'],',');
end

% Generate mound centroid, in case they don't exist already
[geom, ~] = polygeom(m(container).m(mound).X(1:end-2), m(container).m(mound).Y(1:end-2));
mat(container).mcent(mound).mc = [geom(2) , geom(3)];
    
%for nmic = 2:10
%    if length(m(container).m)>=nmic
%        [geom, centr, stuff] = polygeom(m(container).m(nmic).X(1:end-2), m(container).m(nmic).Y(1:end-2))
%        mat(container).mcent(nmic).mc = [geom(2) , geom(3)];
%    end
%end

% Divide the shapefile into separate bed traces by Id
for i = 1:length(newShape)
    divS(newShape(i).Id + 1).mx = [divS(newShape(i).Id + 1).mx newShape(i).X];
    divS(newShape(i).Id + 1).my = [divS(newShape(i).Id + 1).my newShape(i).Y];
    divS(newShape(i).Id + 1).ras = [divS(newShape(i).Id + 1).ras newShape(i).RASTERVALU];
end

for i = 1:length(divS)
    
        shapexval = divS(i).mx;
    shapeyval = divS(i).my;
    shaperval = divS(i).ras;
    
% Define the regress plane for each trace/Id
    [slope(i), azimuth(i), error(i), stdr(i), BVec(i).B, sllb(i), slub(i), azlb(i), azub(i)] =...
        regressplanemod(divS(i).mx',divS(i).my',divS(i).ras');

% Find average x and y coordiantes for each trace/Id
    xmeters(i) = mean(divS(i).mx);
    ymeters(i) = mean(divS(i).my);
    
% Fit calculated xy coords to existing polar coordinates
    dist(i) = real(sqrt((stretchx*(xmeters(i)-mat(container).mcent(mound).mc(1)))^2+(ymeters(i)-mat(container).mcent(mound).mc(2))^2));
    angle(i) = radtodeg(atan2(sqrt((ymeters(i)-mat(container).mcent(mound).mc(2))^2),sqrt((xmeters(i)-mat(container).mcent(mound).mc(1))^2)));
    
%Find the distance to the closest container edge
    contX = c(container).c.X(~isnan(c(container).c.X));
    contY = c(container).c.Y(~isnan(c(container).c.Y));
    v1 = [contX(1:end-1) ; contY(1:end-1)];
    v2 = [v1(1:2,end), v1(1:2,1:end-1)];
    for l = 1:length(v1)
        distc(i) = min(distc(i),distancefunctions([xmeters(i) ymeters(i)]',v1(:,l),v2(:,l), stretchx));
    end
    
% Find the interpolated mound thickness at the traced coordinates
    xvect = linspace(r(container).r.XWorldLimits(1),r(container).r.XWorldLimits(2),length(a(container).a(1,:)));
    yvect = linspace(r(container).r.YWorldLimits(2),r(container).r.YWorldLimits(1),length(a(container).a(:,1)));
    [XV,YV] = meshgrid(xvect,yvect);
    thick(i) = interp2(XV,YV,mat(container).diffM(mound).d,xmeters(i),ymeters(i),'linear');
    elevation(i) = interp2(XV,YV,a(container).a,xmeters(i),ymeters(i),'linear');

% Append any new data points to 'polarDataC[c]M[m].csv'
    match = false;
    if ~isempty(allFiles)
        for j = 1:length(allFiles(:,1))
            if (abs(dist(i) - allFiles(j,3)) <= 1) && (abs(angle(i) - allFiles(j,4)) <= 1E-4)
                match = true
            end
        end
    end
    if match == false
        allFiles = [allFiles;xmeters(i),ymeters(i),dist(i),angle(i),distc(i),thick(i),elevation(i),slope(i),azimuth(i),...
            sllb(i),slub(i),azlb(i),azub(i)];
    end
end
csvwrite(['polarData_C',num2str(container),'M',num2str(mound),'.csv'],allFiles);

figure
contour(r(container).r.XLimWorld(1)+[1:r(container).r.XLimIntrinsic(2)]*r(container).r.DeltaX',r(container).r.YLimWorld(2)+...
    [1:r(container).r.YLimIntrinsic(2)]*r(container).r.DeltaY',a(container).a)
hold on
scatter(allFiles(:,1),allFiles(:,2),100,'blue','filled')
scatter(mat(container).mcent(mound).mc(1),mat(container).mcent(mound).mc(2),100,'red','filled')

figure
errorbar(azimuth,slope,abs(slope-sllb),abs(slope-slub),'x')
hold on
herrorbar(azimuth,slope,abs(azimuth-azlb),abs(azimuth-azub),'x')
for i = 1:length(azimuth)
    text(azimuth(i) + .1*azimuth(i),slope(i) + .1*slope(i),int2str(i))
end