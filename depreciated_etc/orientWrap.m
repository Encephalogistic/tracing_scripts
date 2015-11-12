function [slope, azimuth, error, standErr, BVec] = orientWrap(fileName)

mx = [];
my = [];
ras = [];
divS = struct('mx',mx,'my',my,'ras',ras);

newShape = shaperead(fileName);

for i = 1:length(newShape)
    divS(newShape.Id + 1).mx = [divS(newSHape.Id + 1).mx newShape(i).X];
    divS(newShape.Id + 1).my = [divS(newSHape.Id + 1).my newShape(i).Y];
    divS(newShape.Id + 1).ras = [divS(newSHape.Id + 1).ras newShape(i).RASTERVALU];
end

for i = 1:length(divS)
    [slope(i), azimuth(i), error(i), standErr(i), BVec(i)] = regressplane(divS(i).mx',divS(i).my',divS(i).ras');
end