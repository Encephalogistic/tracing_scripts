ccent = struct('cc',{[]});
mcent = struct('mc',{[]});
mpolarr = struct('mp',{[]});
mpolara = struct('mp',{[]});

%% Calculate the centroid of each container and mound.  Stolen without
%  remorse from HJ Sommer III (polygeom.m, see documentation)

for i = 1:numfiles
    ccent(i).cc = [];
    mcent(i).mc = [];
    [geom, centr, stuff] = polygeom(c(i).c.X(1:end-2), c(i).c.Y(1:end-2))
    ccent(i).cc = [geom(2), geom(3)];
    [geom, centr, stuff] = polygeom(m(i).m.X(1:end-2), m(i).m.Y(1:end-2))
    mcent(i).mc = [geom(2) , geom(3)];
end

%% Calculate mound polar coordinates given relative z.
%  'diffM' is already a measure of z', so it's just a matter of 
%  assigning each pixel a radius and distance from the centroid as well.
%  Assignment is [distance, angle, z']

for i = 1:numfiles
    gpsM = zeros(r(i).r.RasterSize);
    [xs, ys] = size(gpsM);
    gpsM = repmat(gpsM,1,1,3);                 
    stretchx = cos(mcent(i).mc(2)/3376200.00386);
    dx = [r(i).r.XWorldLimits(1):r(i).r.CellExtentInWorldX:r(i).r.XWorldLimits(2)];
    dy = [r(i).r.YWorldLimits(1):r(i).r.CellExtentInWorldY:r(i).r.YWorldLimits(2)];
    for j = 1:length(dx)-1
        for k = 1:length(dy)-1
            %gpsM(j,k,1) = sqrt((stretchx*(dx(j)-mcent(i).mc(1)))^2+(dy(k)-mcent(i).mc(2))^2);
            gpsM(j,k,1) = sqrt( (stretchx*(dx(j)-mcent(i).mc(1)))^2 + (dy(k)-mcent(i).mc(2))^2 );
            gpsM(j,k,2) = atan2(sqrt((dy(k)-mcent(i).mc(2))^2), sqrt((dx(j)-mcent(i).mc(1))^2));
            gpsM(j,k,3) = diffM(i).d(k,j);
        end
    end
    mpolarr(i).mp = gpsM
end

%% Calculate mound polar coordinates given absolute z.
%  The zero will be unique to each mound, so just take min of the
%  interpolated surface.

for i = 1:numfiles
    gpsM = zeros(r(i).r.RasterSize);
    [xs, ys] = size(gpsM);
    gpsM = repmat(gpsM,1,1,3);                 
    stretchx = cos(mcent(i).mc(2)/3376200.00386);
    dx = [r(i).r.XWorldLimits(1):r(i).r.CellExtentInWorldX:r(i).r.XWorldLimits(2)];
    dy = [r(i).r.YWorldLimits(1):r(i).r.CellExtentInWorldY:r(i).r.YWorldLimits(2)];
    minz = min(min(wideInterp(i).w));
    for j = 1:length(dx)-1
        for k = 1:length(dy)-1
            gpsM(j,k,1) = sqrt((stretchx*(dx(j)-mcent(i).mc(1)))^2+(dy(k)-mcent(i).mc(2)));
            gpsM(j,k,2) = atan(sqrt((dy(k)-mcent(i).mc(2))^2)/sqrt((dx(j)-mcent(i).mc(1))^2));
            gpsM(j,k,3) = wideInterp(i).w(k,j)-minz;
        end
    end
    mpolara(i).mp = gpsM
end

% 3376200.00386 meters per radian (or: horizontal distance * cos(YWorldLimits(1)/ <- ))