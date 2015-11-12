function [] = rot_orientation(filename,container,mound)
%============
% File inputs:
%   fileName    => the name of the shapefile group containing layer trace
%                   (required)
%   container   => 1-40, the container surrounding the trace
%   mound       => 1-10, identifies which mound in the container to use as
%                   origin
%============
% Products:
% A contour map of bedding planes within their geologic setting
% A scatter plot of planes graphing slope against azimuth
% A csv in the native directory, named 'bedplane_...', with verbose data
%
%============

load ('MOUND_STATS_WORKING.mat','mat','ccent','c','a','r','m')
[pathname,shapename,~] = fileparts(filename);
traceShape = shaperead(filename);

%Defaults to Mt. Sharp in Gale Crater
if container == 0
    container = 8;
end

if mound == 0
    mound = 1;
end

%Before declaring variables, find out max Id value for preallocation
maxId = 0;
for i = 1:length(traceShape)
    maxId = max(maxId,traceShape(i).ORIG_FID)
end

%maxId = traceShape(length(traceShape)).Id

%Initialize variables for preallocation & convenience
mx = [];
my = [];
ras = [];
divS = struct('mx',mx,'my',my,'ras',ras);

divS(maxId+1).mx = [];  %Divides the shapefile in to separate vectors by bed (equal-length xyz vectors, one for each point in the bed trace)
slope(maxId+1) = 0;     %Slope and azimuth of best-fit plane
azimuth(maxId+1) = 0;
sllb(maxId+1) = 0;
slub(maxId+1) = 0;
azlb(maxId+1) = 0;
azub(maxId+1) = 0;
error(maxId+1) = 0;     %Error measured (currently) by the 'cone of error' method
stdr(maxId+1) = 0;
BVec(maxId+1).B = [];   %The coefficients of each best-fit plane equation
xmeters(maxId+1) = 0;   %XY position of the individual beds (average of the points)
ymeters(maxId+1) = 0;
angle(maxId+1) = 0;     %Radial position relative to mound centroid
dist(maxId+1) = 0;      
thick(maxId+1) = 0;     %Thickness of the mound above interpolated basal surface
elev(maxId+1) = 0;      %Stratigraphic elevation interpolated from MOLA data
radius(maxId+1) = 0;

stretchx = cos(ccent(container).cc(2)/3376200.00386);

for i=1:maxId+1
    distc(i) = nan;     %Distance from nearest container wall
end

%read in old data, or create a csv if a proper one is not found
if exist([pathname,'\bedplane_C',num2str(container),'M',num2str(mound),'.csv'],'file')
    orientationData = csvread([pathname,'\bedplane_C',num2str(container),'M',num2str(mound),'.csv'],1,0);
else
    fid = fopen([pathname,'\bedplane_C',num2str(container),'M',num2str(mound),'.csv'],'w');
    fprintf(fid,'%s\r\n',...
        ['ID,X Coord,Y Coord,Centroid Radius,Centroid Angle,Wall Distance,Mound Thickness,Elevation,',...
        'Slope,UB,LB,Azimuth,UB,LB']);
    fclose(fid);
    ['File created: ',pathname,'\bedplane_C',num2str(container),'M',num2str(mound),'.csv']
    orientationData = [];
end
    
% Divide the shapefile into separate bed traces by Id
for i = 1:length(traceShape)
    divS(traceShape(i).ORIG_FID + 1).mx = [divS(traceShape(i).ORIG_FID + 1).mx traceShape(i).X];
    divS(traceShape(i).ORIG_FID + 1).my = [divS(traceShape(i).ORIG_FID + 1).my traceShape(i).Y];
    divS(traceShape(i).ORIG_FID + 1).ras = [divS(traceShape(i).ORIG_FID + 1).ras traceShape(i).RASTERVALU];
end

% Generate mound centroid
[geom, ~] = polygeom(m(container).m(mound).X(1:end-2), m(container).m(mound).Y(1:end-2));
mat(container).mcent(mound).mc = [geom(2) , geom(3)];

for i = 1:length(divS)
%for i = 1
    
% Define the regress plane for each trace/Id
    [slope(i), azimuth(i), error(i), stdr(i), BVec(i).B, sllb(i), slub(i), azlb(i), azub(i), BINT, R, RINT, STATS] =...
        regressplanemod_v2(divS(i).mx',divS(i).my',divS(i).ras');
    [slope(i), azimuth(i), error(i), stdr(i), BVec(i).B, sllb(i), slub(i), azlb(i), azub(i),~] =...
        regressplanecone(divS(i).mx',divS(i).my',divS(i).ras');
    
    %EDWIN EDIT: Find the point-by-point elevation residuals for the best-fit plane.
    for permutationindex = 1:length(divS(i).mx)
        permutedz = BVec(i).B(3) + BVec(i).B(1).*divS(i).mx' + BVec(i).B(2).*divS(i).my' + circshift(R,permutationindex);
        [slopep(i,permutationindex), azimuthp(i,permutationindex), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] =...
            regressplanemod_v2(divS(i).mx',divS(i).my',permutedz);       
    end

% Find average x and y coordiantes for each trace/Id
    xmeters(i) = mean(divS(i).mx);
    ymeters(i) = mean(divS(i).my);
    
% Fit calculated xy coords to calculated polar coordinates
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
    elev(i) = interp2(XV,YV,a(container).a,xmeters(i),ymeters(i),'linear');

% Append any new data points to 'bedplane_C[c]M[m].csv'
    if isempty(orientationData) || isempty(orientationData(orientationData(:,1) == i-1))
        newData = [double(i-1),xmeters(i),ymeters(i),dist(i),angle(i),distc(i),thick(i),elev(i),...
            slope(i),slub(i),sllb(i),azimuth(i),azub(i),azlb(i)];
        dlmwrite([pathname,'\bedplane_C',num2str(container),'M',num2str(mound),'.csv'], newData,'-append','delimiter',',','roffset',0,'coffset',0);
        orientationData = [orientationData;newData];
    end
end

%Output a contour map with beds located on it
%figure
%contour(r(container).r.XLimWorld(1)+[1:r(container).r.XLimIntrinsic(2)]*r(container).r.DeltaX',r(container).r.YLimWorld(2)+...
%    [1:r(container).r.YLimIntrinsic(2)]*r(container).r.DeltaY',a(container).a)
%hold on
%scatter(orientationData(:,2),orientationData(:,3),100,'blue','filled')
%scatter(mat(container).mcent(mound).mc(1),mat(container).mcent(mound).mc(2),100,'red','filled')

%Plot slope and azimuth with error bars 
figure
errorbar(azimuth,slope,abs(slope-sllb),abs(slope-slub),'x')
hold on
herrorbar(azimuth,slope,abs(azimuth-azlb),abs(azimuth-azub),'x')

slopep(slopep==0)=NaN;azimuthp(azimuthp==0)=NaN;
errorbar(azimuth,slope,slope-min(slopep'),max(slopep')-slope,'x','Color','k')
hold on
h = herrorbar(azimuth,slope,azimuth-min(azimuthp'),max(azimuthp')-azimuth)
set(h,'Color','k')
set(h(2),'LineStyle','none')

for i = 1:length(azimuth)
    text(azimuth(i) + .1*azimuth(i),slope(i) + .1*slope(i),int2str(i-1))
end
