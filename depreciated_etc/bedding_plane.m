%function [relev, angle, dist, slope, azimuth, xmeters, ymeters] = traceToPolar(fileName,container,mound,mat,r,a)

%============
% File inputs:
%   fileName    => the name of the shapefile group containing layer trace
%============

fileName = 'GCLT_1750_VTP';

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
    divS(newShape(i).Id + 1).ras = [divS(newShape(i).Id + 1).ras newShape(i).RASTERVALU];
end

%read in old data
if exist(['polarData_C',num2str(container),'M',num2str(mound),'.csv'],'file')
    allFiles = importdata(['polarData_C',num2str(container),'M',num2str(mound),'.csv'],',');
else
    csvwrite(['polarData_C',num2str(container),'M',num2str(mound),'.csv'],[]);
    ['File created: polarData_C',num2str(container),'M',num2str(mound),'.csv']
    allFiles = importdata(['polarData_C',num2str(container),'M',num2str(mound),mound,'.csv'],',');
end

for i = 1:length(divS)

% Define the regress plane for each trace/Id
    [slope(i), azimuth(i), error(i), stdr(i), BVec(i).B, sllb(i), slub(i), azlb(i), azub(i)] =...
        regressplanemod(divS(i).mx',divS(i).my',divS(i).ras');

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
csvwrite(['bedding_plane_data_',filename,'.csv'],allFiles);

%Plot slope and azimuth with error bars 
figure
errorbar(azimuth,slope,abs(slope-sllb),abs(slope-slub),'x')
hold on
herrorbar(azimuth,slope,abs(azimuth-azlb),abs(azimuth-azub),'x')
for i = 1:length(azimuth)
    text(azimuth(i) + .1*azimuth(i),slope(i) + .1*slope(i),int2str(i))
end