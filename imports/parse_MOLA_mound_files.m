%parse_MOLA_mound_files

numfiles = 40; 

%Ingest files
for l = 1:numfiles
l
name3 = ['container_' int2str(l) '_mola128.tif'];
name3
[a(l).a,r(l).r] = geotiffread(name3);
a(l).a(a(l).a<-2e4) = NaN;
c(l).c = shaperead(['container_' num2str(l)])
m(l).m = shaperead(['mounds_container_' num2str(l)])
end


%%Display mounds
for  nm = 1:numfiles
    figure
contour(r(nm).r.XLimWorld(1)+[1:r(nm).r.XLimIntrinsic(2)]*r(nm).r.DeltaX',r(nm).r.YLimWorld(2)+ ...
    [1:r(nm).r.YLimIntrinsic(2)]*r(nm).r.DeltaY',a(nm).a)
hold on
line(c(nm).c(1).X,c(nm).c(1).Y,'Color','m','LineWidth',2)
line(m(nm).m(1).X,m(nm).m(1).Y,'Color','k','LineWidth',2)

for nmic = 2:10 %number of mounds in container
if size(m(nm).m,1)>=nmic
    line(m(nm).m(nmic).X,m(nm).m(nmic).Y,'Color','k','LineWidth',2)
end
end

end

save MOLA_mounds_working