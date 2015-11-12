%v1 by Jonathan Sneed, Fall 2015
%v2 with revisions by Edwin Kite, October 2015.

%function [relev, angle, dist, slope, azimuth, xmeters, ymeters] = traceToPolar(fileName,container,mound,mat,r,a)

%============
% File inputs:
%   fileName    => the name of the shapefile group containing layer trace
%============

fileName = 'PSP_006855_1750_PSP_007501_1750_trace_VTP';

%============

%Declarations and preallocation
mx = [];
my = [];
ras = [];
maxId = 0;
divS(maxId+1).mx = [];
slope = [];
azimuth = [];
error = [];
stdr = [];
sllb = [];
slub = [];
azlb = [];
azub = [];
divS = struct('mx',mx,'my',my,'ras',ras);

newShape = shaperead(fileName);

for i = 1:length(newShape)
    maxId = max(maxId,newShape(i).Id);
end

divS(maxId+1).mx = [];
slope(maxId+1) = 0;
azimuth(maxId+1) = 0;
error(maxId+1) = 0;
stdr(maxId+1) = 0;
BVec(maxId+1).B = 0;
sllb(maxId+1) = 0;
slub(maxId+1) = 0;
azlb(maxId+1) = 0;
azub(maxId+1) = 0;

% Divide the shapefile into separate bed traces by Id
for i = 1:length(newShape)
    divS(newShape(i).Id + 1).mx = [divS(newShape(i).Id + 1).mx newShape(i).X];
    divS(newShape(i).Id + 1).my = [divS(newShape(i).Id + 1).my newShape(i).Y];
    divS(newShape(i).Id + 1).ras = [divS(newShape(i).Id + 1).ras newShape(i).RASTERVALU]; %elevation values
end

%read in old data
if exist(['bedding_plane_data_',fileName,'.csv'],'file')
    allFiles = importdata(['bedding_plane_data_',fileName,'.csv'],',');
else
    csvwrite(['bedding_plane_data_',fileName,'.csv'],[]);
    ['File created: bedding_plane_data_',fileName,'.csv']
    allFiles = importdata(['bedding_plane_data_',fileName,'.csv'],',');
end

for i = 1:length(divS) %For each bed trace

% Define the regress plane for each trace/Id
    [slope(i), azimuth(i), error(i), stdr(i), BVec(i).B, sllb(i), slub(i), azlb(i), azub(i), BINT, R, RINT, STATS] =...
        regressplanemod_v2(divS(i).mx',divS(i).my',divS(i).ras');

%EDWIN EDIT: Find the point-by-point elevation residuals for the best-fit plane.
    for permutationindex = 1:length(divS(i).mx)
       permutedz = BVec(i).B(3) + BVec(i).B(1).*divS(i).mx' + BVec(i).B(2).*divS(i).my' + circshift(R,permutationindex);
    [slopep(i,permutationindex), azimuthp(i,permutationindex), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] =...
        regressplanemod_v2(divS(i).mx',divS(i).my',permutedz);       
    end
       
    
% Append any new data points to 'polarDataC[c]M[m].csv'
    match = false;
    if ~isempty(allFiles)
        for j = 1:length(allFiles(:,1))
            if (abs(slope(i) - allFiles(j,3)) <= 1E-4) && (abs(angle(i) - allFiles(j,4)) <= 1E-4)
                match = true
            end
        end
    end
    if match == false
        allFiles = [allFiles;slope(i),sllb(i),slub(i),azimuth(i),azlb(i),azub(i)];
    end
end
csvwrite(['bedding_plane_data_',fileName,'.csv'],allFiles);

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
    text(azimuth(i) + .1*azimuth(i),slope(i) + .1*slope(i),int2str(i))
end