%% Major variable declarations

discX = struct('d',{[]},'onp',{[]});    % ordered x-coordinates for the mound polygon vertices and edges
discY = struct('d',{[]},'onp',{[]});    % ordered y-coordinates for the mound polygon vertices and edges
moundA = struct('m',{[]});              % new matrix generated zeroing the exterior of a mound shape
maxX = 0;                               % smaller version of moundA for efficiency reasons
maxY = 0;                               % contains positioning data for placement within moundA
minX = 0;
minY = 0;
moundAred = struct('m',{[]}, 'maxx', maxX, 'maxy', maxY, 'minx', minX, 'miny', minY);
index = struct('i',{[]});               % Stores selected a values for interpolating 
surfacecube = struct ('s',{[]});        % The interpolated surface, in a's coordinates
wideInterp = struct('w',{[]}); 
diffM = struct('d',{[]});
moundVolume = struct('m',0);
closest = struct('c',{[]});             %Stores the shortest distance to any container wall for each mound point.
mdistc = struct('md',{[]});             %Mound area within a given distance
cudist = struct('c',{[]});              %(Distance is implicit in vector index)
cudist2 = struct('c2',{[]});            %Mound area / container area, by distance
moundWidth = struct('mw',0);
moundRelief = struct('mr',0);
mcent = struct('mc',{[]});
mpolarr = struct('mp',{[]});
mpolara = struct('mp',{[]});

% 'Master' mound attribute struct, mat(M).moundA(N).m is the Nth mound in
% the Mth container, and so on.
mat = struct('discX',{discX},'discY',{discY},'moundA',{moundA},'moundAred',{moundAred},...
    'wideInterp',{wideInterp},'diffM',{diffM},'moundVolume',{moundVolume},...
    'index',{index},'surfacecube',{surfacecube},'closest',{closest},'mdistc',{mdistc},...
    'mcent',{mcent},'mpolara',{mpolara},'mpolarr',{mpolarr},...
    'cudist',{cudist},'cudist2',{cudist2},'moundWidth',{moundWidth},'moundRelief',{moundRelief},...
    'moundAF',{[]},'closestF',{[]},'cudistF',{[]},'cudist2F',{[]},'moundWidthF',0,'moundRF',0);

containerArea = struct('ca',{[]});
containerRelief = struct('cr',{[]});
cdistc = struct('cd',{[]});             %Container area within a given distance
cclosest = struct('cc',{[]});           %Stores the shortest distance to any container wall for each container point.
gps = struct('g',{[]});                 %Links absolute matrix to raster values- gives GPS coordinates of each cell
                                        %(i.e. gps(i).g(5,10) gives X and Y
                                        % for a(i).a(5,10))
    

numfiles = 8;                          %Change from 1-40, leave at 40 if not testing

%% Ingest files

for l = 1:numfiles
    l
    name3 = ['container_' int2str(l) '_mola128.tif'];
    name3
    [a(l).a,r(l).r] = geotiffread(name3);
    a(l).a(a(l).a<-2e4) = NaN;
    c(l).c = shaperead(['container_' num2str(l)])
    m(l).m = shaperead(['mounds_container_' num2str(l)])
end


%% Display mounds

for  nm = 1:numfiles
    figure
    contour(r(nm).r.XLimWorld(1)+[1:r(nm).r.XLimIntrinsic(2)]*r(nm).r.DeltaX',r(nm).r.YLimWorld(2)+[1:r(nm).r.YLimIntrinsic(2)]*r(nm).r.DeltaY',a(nm).a)
    hold on
    line(c(nm).c(1).X,c(nm).c(1).Y,'Color','m','LineWidth',2)
    line(m(nm).m(1).X,m(nm).m(1).Y,'Color','k','LineWidth',2)

    for nmic = 2:10 %number of mounds in container
        if size(m(nm).m,1)>=nmic
            line(m(nm).m(nmic).X,m(nm).m(nmic).Y,'Color','k','LineWidth',2)
        end
    end
    title(['Container ' num2str(nm)])
    xlabel('Longitude')
    ylabel('Latitude')
end

%% Recast shapefile info in discrete coordinates, and find 
% polygon maximums/minimums appropriately.  (Ordered pairs discX, discY
% list the vertices of each mound outline shape.)

for i = 1:numfiles
    mat(i).discX(1).d = zeros(size(m(i).m(1).X));
    mat(i).discY(1).d = zeros(size(m(i).m(1).Y));
   [mat(i).discX(1).d,mat(i).discY(1).d] = worldToDiscrete(r(i).r,m(i).m(1).X,m(i).m(1).Y); 
   mat(i).moundAred(1).maxx = max(mat(i).discX(1).d);
   mat(i).moundAred(1).maxy = max(mat(i).discY(1).d);
   mat(i).moundAred(1).minx = min(mat(i).discX(1).d);
   mat(i).moundAred(1).miny = min(mat(i).discY(1).d);
   
   for nmic = 2:10
        if length(m(i).m)>=nmic
            mat(i).discX(nmic).d = zeros(size(m(i).m(nmic).X));
            mat(i).discY(nmic).d = zeros(size(m(i).m(nmic).Y));
            [mat(i).discX(nmic).d,mat(i).discY(nmic).d] =...
                worldToDiscrete(r(i).r,m(i).m(nmic).X,m(i).m(nmic).Y); 
            mat(i).moundAred(nmic).maxx = max(mat(i).discX(nmic).d);
            mat(i).moundAred(nmic).maxy = max(mat(i).discY(nmic).d);
            mat(i).moundAred(nmic).minx = min(mat(i).discX(nmic).d);
            mat(i).moundAred(nmic).miny = min(mat(i).discY(nmic).d);
        end
    end
end

%% Generate the moundA and moundAred matrices as zeros, then populate it with values only
% from within the given mound shape, leaving everything else as 0.
% Populates the ordered list of points on the polygon vertices for later,
% nefarious uses.  (i.e. surface_gen files)

for i = 1:numfiles
   i
   mat(i).moundA(1).m = zeros(size(a(i).a));
   mat(i).moundAred(1).m = zeros(mat(i).moundAred(1).maxx+2 ...
        -mat(i).moundAred(1).minx,mat(i).moundAred(1).maxy+2-mat(i).moundAred(1).miny);
   [e,f]=size(mat(i).moundAred(1).m);
   
   for j = 1:e
       for k = 1:f
           [inp, onp] = inpolygon(j-1+mat(i).moundAred(1).minx,k-1+...
               mat(i).moundAred(1).miny,mat(i).discX(1).d,mat(i).discY(1).d);
           if inp || onp
             mat(i).moundAred(1).m(j,k) = a(i).a(j-1+mat(i).moundAred(1).minx,k-1+mat(i).moundAred(1).miny);
             mat(i).moundA(1).m(j-1+mat(i).moundAred(1).minx,k-1+mat(i).moundAred(1).miny) =...
                 a(i).a(j-1+mat(i).moundAred(1).minx,k-1+mat(i).moundAred(1).miny);
           end
       end
   end
   mat(i).moundAF = mat(i).moundA(1).m;
   
   for nmic = 2:10
        if length(m(i).m)>=nmic
            mat(i).moundA(nmic).m = zeros(size(a(i).a));
            mat(i).moundAred(nmic).m = zeros(mat(i).moundAred(nmic).maxx+2 ...
                -mat(i).moundAred(nmic).minx,mat(i).moundAred(nmic).maxy+2-mat(i).moundAred(nmic).miny);
            [e,f]=size(mat(i).moundAred(nmic).m);
            
            for j = 1:e
                for k = 1:f
                    [inp, onp] = inpolygon(j-1+mat(i).moundAred(nmic).minx,k-1+...
                        mat(i).moundAred(nmic).miny,mat(i).discX(nmic).d,mat(i).discY(nmic).d);
                    if inp || onp
                        mat(i).moundAred(nmic).m(j,k) = a(i).a(j-1+mat(i).moundAred(nmic).minx,k-1+mat(i).moundAred(nmic).miny);
                        mat(i).moundA(nmic).m(j-1+mat(i).moundAred(nmic).minx,k-1+mat(i).moundAred(nmic).miny) =...
                            a(i).a(j-1+mat(i).moundAred(nmic).minx,k-1+mat(i).moundAred(nmic).miny);
                    end
                end
            end
            mat(i).moundAF = mat(i).moundAF + mat(i).moundA(nmic).m;
        end
    end
end

%% Sample points.
%  Gather points to be used for anchoring the inerpolated surface.  At the
%  moment, this is just the vertices themselves.

for j = 1:numfiles
   mat(j).index(1).i = [];
   for k = 1:length(mat(j).discX(1).d(1:end-1))
       mat(j).index(1).i = [mat(j).index(1).i a(j).a( mat(j).discX(1).d(k), mat(j).discY(1).d(k))];
   end
   
   for nmic = 2:10
        if length(m(j).m)>=nmic
            mat(j).index(nmic).i = [];
            for k = 1:length(mat(j).discX(nmic).d(1:end-1))
                mat(j).index(nmic).i = [mat(j).index(nmic).i a(j).a( mat(j).discX(nmic).d(k), mat(j).discY(nmic).d(k))];
            end
        end
    end
end

%% Interpolation
%  Generate a smooth basal surface from points around the edge of the mound
%  polygon.

for j = 1:numfiles
    [xsize,ysize] = size(a(j).a);
    
    mat(j).surfacecube(1).s = griddata(mat(j).discX(1).d(1:end-1), mat(j).discY(1).d(1:end-1),...
        double(mat(j).index(1).i),1:xsize,(1:ysize)', 'cubic');
    mat(j).surfacecube(1).s(isnan(mat(j).surfacecube(1).s)) = 0; 
    
    for nmic = 2:10
        if length(m(j).m)>=nmic
            mat(j).surfacecube(nmic).s = griddata(mat(j).discX(nmic).d(1:end-1), mat(j).discY(nmic).d(1:end-1),...
                double(mat(j).index(nmic).i),1:xsize,(1:ysize)', 'cubic');
            mat(j).surfacecube(nmic).s(isnan(mat(j).surfacecube(nmic).s)) = 0; 
        end
    end
end

%% Create a container surface with a hole where the mound used to be.
%  Then, fill it with the interpolated surface.  Basically, 'subtract the
%  mound' from the container.

for i = 1:numfiles
    inverseA = double(mat(i).moundA(1).m == 0) .* double(a(i).a);
    mat(i).wideInterp(1).w = inverseA;
    inZone = inverseA == 0;
    mat(i).wideInterp(1).w = mat(i).wideInterp(1).w + inZone .* mat(i).surfacecube(1).s';
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            inverseA = double(mat(i).moundA(nmic).m == 0) .* double(a(i).a);
            mat(i).wideInterp(nmic).w = inverseA;
            inZone = inverseA == 0;
            mat(i).wideInterp(nmic).w = mat(i).wideInterp(nmic).w + inZone .* mat(i).surfacecube(nmic).s';
        end
    end
end

%% Create a new matrix of the difference between observed subaerial surface
%  and interpolated basal surface at each matrix coordinate.  Some of these
%  distances will be negative, where the interpolated surface has a higher
%  elevation than the observed surface of the mound.

for i = 1:numfiles
    mat(i).diffM(1).d = double(a(i).a) - double(mat(i).wideInterp(1).w);
    mat(i).diffM(1).d(mat(i).diffM(1).d == 0) = NaN;
    mat(i).moundVolume(1).m = sum(sum(mat(i).diffM(1).d));
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            mat(i).diffM(nmic).d = double(a(i).a) - double(mat(i).wideInterp(nmic).w);
            mat(i).diffM(nmic).d(mat(i).diffM(nmic).d == 0) = NaN;
            mat(i).moundVolume(nmic).m = sum(sum(mat(i).diffM(nmic).d));
        end
    end
end

%% Calculate the centroid of each container.  Stolen without
%  remorse from HJ Sommer III (polygeom.m, see documentation).  Then,
%  assemble a vector assigning each point in the container and mound to a
%  distance value.  The vector doesn't record the original index in the
%  matrix, but preserves the left-to-right up-to-down reading order for
%  later reconstruction if desired (I do this for the distance plots at the
%  end).

for i = 1:numfiles
    contX = c(i).c.X(~isnan(c(i).c.X));
    contY = c(i).c.Y(~isnan(c(i).c.Y));
    
    %Find the centroid of the container.
    [geom, centr, stuff] = polygeom(contX(1:end-1), contY(1:end-1));
    ccent(i).cc = [geom(2), geom(3)];

    %Generate points to be measured with real distances rather than matrix coordiantes.
    gps(i).g = [];
    [p,q] = meshgrid(linspace(r(i).r.XWorldLimits(1),r(i).r.XWorldLimits(2),r(i).r.RasterSize(2)),...
        linspace(r(i).r.YWorldLimits(2),r(i).r.YWorldLimits(1),r(i).r.RasterSize(1)));
    gps(i).g = [p(:) q(:)];
    
    i
    cclosest(i).cc = NaN(1,length(gps(i).g));
    size(cclosest(i).cc)
    
    % Define line segments
    v1 = [contX(1:end-1) ; contY(1:end-1)];
    v2 = [v1(1:2,end), v1(1:2,1:end-1)];
    
    %Calculate distances
    for l = 1:length(v1)
        a_points = repmat(v1(:,l),1,length(gps(i).g));
        b_points = repmat(v2(:,l),1,length(gps(i).g));
        cclosest(i).cc = min(cclosest(i).cc,distancefunctions(gps(i).g', a_points, b_points, cos(ccent(i).cc(2)/3376200.00386)));
    end
    
    tempMatrix = nan(size(reshape(cclosest(i).cc,size(p))));
    
    for k = 1:length(p(:,1))
        for j = 1:length(p(1,:))
            out = false;
            if a(i).a(k,j)==0&&~out
                inpolygon(j-1+mat(i).moundAred(1).minx,k-1+mat(i).moundAred(1).miny,c(i).c.X,c(i).c.Y);
                out = true;
            else
                out = false;
                tempMatrix = 1;
            end                
        end
    end
    cclosest(i).cc = cclosest(i).cc .* reshape(double(tempMatrix),size(cclosest(i).cc));
    %cclosest(i).cc(cclosest(i).cc == 0) = NaN;
    
    mat(i).closest(1).c = cclosest(i).cc .* reshape(double(mat(i).moundA(1).m~=0),size(cclosest(i).cc));
    mat(i).closest(1).c(mat(i).closest(1).c == 0) = NaN;
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            mat(i).closest(nmic).c = cclosest(i).cc .* reshape(double(mat(i).moundA(nmic).m~=0),size(cclosest(i).cc));
            mat(i).closest(nmic).c(mat(i).closest(nmic).c == 0) = NaN;
        end
    end
    
    mat(i).closestF = cclosest(i).cc .* reshape(double(mat(i).moundAF~=0),size(cclosest(i).cc));
    mat(i).closestF(mat(i).closestF == 0) = NaN;

    figure
    hold on;
    surf(reshape(cclosest(i).cc,size(p)), 'EdgeColor','none');
    title(['Container ' num2str(i)])
end

%% Display naive cdf area plots

%metaVector = [];
%figure
%hold on;


%for i = 1:numfiles
%    cdfplot(reshape(mat(i).closest(1).c,1,[]))
%    metaVector = [metaVector reshape(mat(i).closest(1).c,1,[])];
    
%    for nmic = 2:10
%        if length(m(i).m)>=nmic
%            cdfplot(reshape(mat(i).closest(nmic).c,1,[]))
%            metaVector = [metaVector reshape(mat(i).closest(nmic).c,1,[])];
%        end
%    end
%end

%figure
%cdfplot(metaVector)

%% Construct the various distribution functions.
%  Does so in a fairly blunt way.  Sort every cell by distance, iterate
%  from zero to the maximum distance, and record the number of cells that
%  fall within each iteration.  Record the delta and rolling totals.

far = 0;
bin = 600;

for i = 1:numfiles
    i
    mi = 1;
    ci = 1;
    fi = 1;
    mdsort = sort(mat(i).closest(1).c);
    fsort = sort(mat(i).closestF);
    cdsort = sort(cclosest(i).cc);
    ttemp1 = 0;
    ttemp2 = 0;
    ttempf = 0;
    
    mat(i).mdistc(1).md = zeros(1,ceil(max(max(cclosest(i).cc))/bin));
    cdistc(i).cd = zeros(1,ceil(max(max(cclosest(i).cc))/bin));
    mat(i).cudist2(1).c2 = zeros(1,ceil(max(max(cclosest(i).cc))/bin));
    
    if ceil(max(max(cclosest(i).cc))/bin) > far
        far = ceil(max(max(cclosest(i).cc))/bin);
    end
        
    for rd = 2:ceil(max(max(cclosest(i).cc))/bin)
        temp1 = 0;
        while mdsort(mi)<rd*bin
            temp1 = temp1 + mdsort(mi);
            mi = mi + 1;
        end
        mat(i).mdistc(1).md(rd) = mat(i).mdistc(1).md(rd-1)+temp1;
        ttemp1 = ttemp1 + temp1;
                
        temp2 = 0;
        while cdsort(ci)<rd*bin
            temp2 = temp2 + cdsort(ci);
            ci = ci + 1;
        end
        cdistc(i).cd(rd) = cdistc(i).cd(rd-1)+temp2;
        ttemp2 = ttemp2 + temp2;     
        
        tempf = 0;
        while fsort(fi)<rd*bin
            tempf = tempf + fsort(fi);
            fi = fi + 1;
        end
        ttempf = ttempf + tempf;     
        
        if temp2 == 0
            mat(i).cudist2(1).c2(rd) = mat(i).cudist2(1).c2(rd-1);
            mat(i).cudist2F(rd) = mat(i).cudist2F(rd-1);
        else
            mat(i).cudist2(1).c2(rd) = temp1/temp2;% + cudist2(i).c2(rd-1);
            mat(i).cudist2F(rd) = tempf/temp2;
        end
        mat(i).cudist(1).c(rd) = ttemp1/ttemp2;
        mat(i).cudistF(rd) = ttempf/ttemp2;
    end
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            mi = 1;
            ci = 1;
            mdsort = sort(mat(i).closest(nmic).c);
            ttemp1 = 0;
            ttemp2 = 0;
            
            mat(i).mdistc(nmic).md = zeros(1,ceil(max(max(cclosest(i).cc))/bin));
            
            for rd = 2:ceil(max(max(cclosest(i).cc))/bin)
                temp1 = 0;
                while mdsort(mi)<rd*bin
                    temp1 = temp1 + mdsort(mi);
                    mi = mi + 1;
                end
                mat(i).mdistc(nmic).md(rd) = mat(i).mdistc(nmic).md(rd-1)+temp1;
                ttemp1 = ttemp1 + temp1;
                
                temp2 = 0;
                while cdsort(ci)<rd*bin
                    temp2 = temp2 + cdsort(ci);
                    ci = ci + 1;
                end
                cdistc(i).cd(rd) = cdistc(i).cd(rd-1)+temp2;
                ttemp2 = ttemp2 + temp2;       
        
                if temp2 == 0
                    mat(i).cudist2(nmic).c2(rd) = mat(i).cudist2(nmic).c2(rd-1);
                else
                    mat(i).cudist2(nmic).c2(rd) = temp1/temp2;% + cudist2(i).c2(rd-1);
                end
                mat(i).cudist(nmic).c(rd) = ttemp1/ttemp2;
            end
        end
    end
end

%% Show cumulative plots, and the aggregate plot.
%  Normalizes all mounds against one another- the aggregate plot weights
%  all mounds equally.

figure
hold on
metadist = zeros(1,far);
for i = 1:numfiles
   mat(i).cudist2(1).c2(isnan(mat(i).cudist2(1).c2)) = max(mat(i).cudist2(1).c2);
   mat(i).cudist2(1).c2 = [mat(i).cudist2(1).c2 ./ max(mat(i).cudist2(1).c2), ones(1,length(metadist)-length(mat(i).cudist2(1).c2))];
   plot((1:length(mat(i).cudist2(1).c2)) .* bin,mat(i).cudist2(1).c2)
   metadist = metadist + [mat(i).cudist2(1).c2, zeros(1,length(metadist)-length(mat(i).cudist2(1).c2))];
   
   for nmic = 2:10
        if length(m(i).m)>=nmic
            mat(i).cudist2(nmic).c2(isnan(mat(i).cudist2(nmic).c2)) = max(mat(i).cudist2(nmic).c2);
            mat(i).cudist2(nmic).c2 = [mat(i).cudist2(nmic).c2 ./ max(mat(i).cudist2(nmic).c2), ones(1,length(metadist)-length(mat(i).cudist2(nmic).c2))];
            plot((1:length(mat(i).cudist2(nmic).c2)) .* bin,mat(i).cudist2(nmic).c2)
            metadist = metadist + [mat(i).cudist2(nmic).c2, zeros(1,length(metadist)-length(mat(i).cudist2(nmic).c2))];
        end
    end
end
title(['Fraction of the container within each mound'])
xlabel('Distance from container wall (/600m)')
ylabel('Fraction within a specific distance bin')

figure
plot((1:length(metadist)) .* bin,metadist./numfiles)
title('Aggregate fraction of container in bin, over all mounds & containers')
xlabel('Distance from container wall (/600m)')
ylabel('Fraction within a specific distance bin')

figure
hold on
for i = 1:numfiles
   plot((1:length(mat(i).cudist(1).c)) .* bin,mat(i).cudist(1).c ./ max(mat(i).cudist(1).c))
   
   for nmic = 2:10
        if length(m(i).m)>=nmic
            plot((1:length(mat(i).cudist(nmic).c)) .* bin,mat(i).cudist(nmic).c ./ max(mat(i).cudist(nmic).c))
        end
    end
end
title('Normalized cumulative fraction of container within the mound')
xlabel('Distance from container wall (/600m)')
ylabel('Mound fraction of container up to specific distance')

figure
hold on
for i = 1:numfiles
    plot((1:length(mat(i).cudist(1).c)) .* bin,mat(i).cudist(1).c)
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            plot((1:length(mat(i).cudist(nmic).c)) .* bin,mat(i).cudist(nmic).c)
        end
    end
end
title('Cumulative fraction of container within the mound')
xlabel('Distance from container wall (/600m)')
ylabel('Mound fraction of container up to specific distance')

%% Show cumulative plots, and the aggregate plot.
%  Same as the previous code bloc, only with all mounds mooshed in to one
%  giant mound per container.

figure
hold on
metadist = zeros(1,far);
for i = 1:numfiles
   mat(i).cudist2F(isnan(mat(i).cudist2F)) = max(mat(i).cudist2F);
   mat(i).cudist2F = [mat(i).cudist2F ./ max(mat(i).cudist2F), ones(1,length(metadist)-length(mat(i).cudist2F))];
   plot(mat(i).cudist2(1).c2)
   metadist = metadist + [mat(i).cudist2F, zeros(1,length(metadist)-length(mat(i).cudist2F))];
end
title(['Fraction of the container within each mound (C.M.)'])
xlabel('Distance from container wall (/600m)')
ylabel('Fraction within a specific distance bin')

figure
plot(metadist./numfiles)
title('Aggregate fraction of container in bin, over all mounds & containers (C.M.)')
xlabel('Distance from container wall (/600m)')
ylabel('Fraction within a specific distance bin')

figure
hold on
for i = 1:numfiles
   plot(mat(i).cudistF ./ max(mat(i).cudistF))
end
title('Normalized cumulative fraction of container within the mound (C.M.)')
xlabel('Distance from container wall (/600m)')
ylabel('Mound fraction of container up to specific distance')

figure
hold on
for i = 1:numfiles
    plot(mat(i).cudistF)
end
title('Cumulative fraction of container within the mound (C.M.)')
xlabel('Distance from container wall (/600m)')
ylabel('Mound fraction of container up to specific distance')

%% Generate moundWidth values

for i = 1:numfiles
    moundW = 0;
    index = 1;
    
    while moundW == 0
        if mat(i).cudist2(1).c2(index) >= 0.5
            moundW = index;
        end
        index = index + 1;
    end
    mat(i).moundWidth(1).mw = moundW;
    
    moundW = 0;
    index = 1;
    
    while moundW == 0
        if mat(i).cudist2F(index) >= 0.5
            moundW = index;
        end
        index = index + 1;
    end
    mat(i).moundWidthF = moundW;
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            moundW = 0;
            index = 1;
            
            while moundW == 0
                if mat(i).cudist2(nmic).c2(index) >= 0.5
                    moundW = index;
                end
                index = index + 1;
            end
            mat(i).moundWidth(nmic).mw = moundW;
        end
    end
end

%% Calculate container area

for i = 1:numfiles
    containerArea(i).ca = sum(sum(~isnan(a(i).a./a(i).a)));
end

 %% Calculate container relief
 
for i = 1:numfiles
    %tempInterp = wideInterp(i).w;
    tempInterp = a(i).a;
    tempInterp(tempInterp  == 0) = NaN;
    containerRelief(i).cr = sqrt((double(max(max(tempInterp)))-double(min(min(tempInterp))))^2);
end

%% Calculate mound relief

for i = 1:numfiles
    mat(i).moundRelief(1).mr = max(max(mat(i).diffM(1).d));
    mat(i).moundRF = mat(i).moundRelief(1).mr;
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            mat(i).moundRelief(nmic).mr = max(max(mat(i).diffM(nmic).d));
            mat(i).moundRF = max(mat(i).moundRF,mat(i).moundRelief(nmic).mr);
        end
    end
end

%% Scatter plotes

figure
hold on
for i = 1:numfiles
    if i == 5 ||...
        i == 6 ||...
        i == 7 ||...
        i == 9||...
        i == 10||...
        i == 12||...
        i == 14
        scatter(containerArea(i).ca,mat(i).moundWidth(1).mw,100,'red','filled')
    elseif i == 4||...
            i == 23
        scatter(containerArea(i).ca,mat(i).moundWidth(1).mw,100,'green','filled')
    else
        scatter(containerArea(i).ca,mat(i).moundWidth(1).mw,100,'blue','filled')
    end
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            if i == 5 ||...
                i == 6 ||...
                i == 7 ||...
                i == 9||...
                i == 10||...
                i == 12||...
                i == 14
                scatter(containerArea(i).ca,mat(i).moundWidth(nmic).mw,100,'red','filled')
            elseif i == 4||...
                i == 23
                scatter(containerArea(i).ca,mat(i).moundWidth(nmic).mw,100,'green','filled')
            else
                scatter(containerArea(i).ca,mat(i).moundWidth(nmic).mw,100,'blue','filled')
            end
        end
    end
end
title('Mound width as a function of container area')
xlabel('Number of pixels per container')
ylabel('Mound width /600m')

figure
hold on
for i = 1:numfiles
    if i == 5 ||...
        i == 6 ||...
        i == 7 ||...
        i == 9||...
        i == 10||...
        i == 12||...
        i == 14
        scatter(containerRelief(i).cr,mat(i).moundWidth(1).mw,100,'red','filled')
    elseif i == 4||...
            i == 23
        scatter(containerRelief(i).cr,mat(i).moundWidth(1).mw,100,'green','filled')
    else
        scatter(containerRelief(i).cr,mat(i).moundWidth(1).mw,100,'blue','filled')
    end
  
    for nmic = 2:10
        if length(m(i).m)>=nmic
            if i == 5 ||...
                i == 6 ||...
                i == 7 ||...
                i == 9||...
                i == 10||...
                i == 12||...
                i == 14
                scatter(containerRelief(i).cr,mat(i).moundWidth(nmic).mw,100,'red','filled')
            elseif i == 4||...
                i == 23
                scatter(containerRelief(i).cr,mat(i).moundWidth(nmic).mw,100,'green','filled')
            else
                scatter(containerRelief(i).cr,mat(i).moundWidth(nmic).mw,100,'blue','filled')
            end
        end
    end
end
title('Mound width as a function of container relief')
xlabel('Container relief (m)')
ylabel('Mound width /600m')

figure
hold on
for i = 1:numfiles
    if i == 5 ||...
        i == 6 ||...
        i == 7 ||...
        i == 9||...
        i == 10||...
        i == 12||...
        i == 14
        scatter(mat(i).moundRelief(1).mr,mat(i).moundWidth(1).mw,100,'red','filled')
    elseif i == 4||...
            i == 23
        scatter(mat(i).moundRelief(1).mr,mat(i).moundWidth(1).mw,100,'green','filled')
    else
        scatter(mat(i).moundRelief(1).mr,mat(i).moundWidth(1).mw,100,'blue','filled')
    end
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            if i == 5 ||...
                i == 6 ||...
                i == 7 ||...
                i == 9||...
                i == 10||...
                i == 12||...
                i == 14
                scatter(mat(i).moundRelief(nmic).mr,mat(i).moundWidth(nmic).mw,100,'red','filled')
            elseif i == 4||...
                i == 23
                scatter(mat(i).moundRelief(nmic).mr,mat(i).moundWidth(nmic).mw,100,'green','filled')
            else
                scatter(mat(i).moundRelief(nmic).mr,mat(i).moundWidth(nmic).mw,100,'blue','filled')
            end
        end
    end
end
title('Mound width as a function of mound relief')
xlabel('Mound relief (m)')
ylabel('Mound width /600m')

%% Scatter plots (C.M.)

figure
hold on
for i = 1:numfiles
    if i == 5 ||...
        i == 6 ||...
        i == 7 ||...
        i == 9||...
        i == 10||...
        i == 12||...
        i == 14
        scatter(containerArea(i).ca,mat(i).moundWidthF,100,'red','filled')
    elseif i == 4||...
            i == 23
        scatter(containerArea(i).ca,mat(i).moundWidthF,100,'green','filled')
    else
        scatter(containerArea(i).ca,mat(i).moundWidthF,100,'blue','filled')
    end
end
title('Mound width as a function of container area (C.M.)')
xlabel('Number of pixels per container')
ylabel('Mound width /600m')

figure
hold on
for i = 1:numfiles
    if i == 5 ||...
        i == 6 ||...
        i == 7 ||...
        i == 9||...
        i == 10||...
        i == 12||...
        i == 14
        scatter(containerRelief(i).cr,mat(i).moundWidthF,100,'red','filled')
    elseif i == 4||...
            i == 23
        scatter(containerRelief(i).cr,mat(i).moundWidthF,100,'green','filled')
    else
        scatter(containerRelief(i).cr,mat(i).moundWidthF,100,'blue','filled')
    end
end
title('Mound width as a function of container relief (C.M.)')
xlabel('Container relief (m)')
ylabel('Mound width /600m')

figure
hold on
for i = 1:numfiles
    if i == 5 ||...
        i == 6 ||...
        i == 7 ||...
        i == 9||...
        i == 10||...
        i == 12||...
        i == 14
        scatter(mat(i).moundRF,mat(i).moundWidthF,100,'red','filled')
    elseif i == 4||...
            i == 23
        scatter(mat(i).moundRF,mat(i).moundWidthF,100,'green','filled')
    else
        scatter(mat(i).moundRF,mat(i).moundWidthF,100,'blue','filled')
    end
 end
title('Mound width as a function of mound relief (C.M.)')
xlabel('Mound relief (m)')
ylabel('Mound width /600m')