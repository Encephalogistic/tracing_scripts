%% Calculate the centroid of each container and mound.  Stolen without
%  remorse from HJ Sommer III (polygeom.m, see documentation)

for i = 1:numfiles
    [geom, centr, stuff] = polygeom(m(i).m(1).X(1:end-2), m(i).m(1).Y(1:end-2))
    mat(i).mcent(1).mc = [geom(2) , geom(3)];
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            [geom, centr, stuff] = polygeom(m(i).m(nmic).X(1:end-2), m(i).m(nmic).Y(1:end-2))
            mat(i).mcent(nmic).mc = [geom(2) , geom(3)];
        end
    end
end

%% Calculate mound polar coordinates given relative z.
%  'diffM' is already a measure of z', so it's just a matter of 
%  assigning each pixel a radius and distance from the centroid as well.
%  Assignment is [distance, angle, z']

for i = 1:numfiles
    gpsM = zeros(r(i).r.RasterSize);
    [xs, ys] = size(gpsM);
    gpsM = repmat(gpsM,1,1,3);                 
    stretchx = cos(ccent(i).cc(2)/3376200.00386);
    dx = [r(i).r.XWorldLimits(1):r(i).r.CellExtentInWorldX:r(i).r.XWorldLimits(2)];
    dy = [r(i).r.YWorldLimits(1):r(i).r.CellExtentInWorldY:r(i).r.YWorldLimits(2)];
    for j = 1:length(dx)-1
        for k = 1:length(dy)-1
            gpsM(length(dy)-k,length(dx)-j,1) = real(sqrt( (stretchx*(dx(j)-mat(i).mcent(1).mc(1)))^2 + (dy(k)-mat(i).mcent(1).mc(2))^2 ));
            gpsM(length(dy)-k,length(dx)-j,2) = atan2(sqrt((dy(k)-mat(i).mcent(1).mc(2))^2), sqrt((dx(j)-mat(i).mcent(1).mc(1))^2));
            gpsM(k,j,3) = mat(i).diffM(1).d(k,j);
        end
    end
    mat(i).mpolarr(1).mp = gpsM;
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            for j = 1:length(dx)-1
                for k = 1:length(dy)-1
                    gpsM(length(dy)-k,length(dx)-j,1) = real(sqrt( (stretchx*(dx(j)-mat(i).mcent(nmic).mc(1)))^2 + (dy(k)-mat(i).mcent(nmic).mc(2))^2 ));
                    gpsM(length(dy)-k,length(dx)-j,2) = atan2(sqrt((dy(k)-mat(i).mcent(nmic).mc(2))^2), sqrt((dx(j)-mat(i).mcent(nmic).mc(1))^2));
                    gpsM(k,j,3) = mat(i).diffM(nmic).d(k,j);
                end
            end
            mat(i).mpolarr(nmic).mp = gpsM;
        end
    end
end

%% Calculate mound polar coordinates given absolute z.
%  The zero will be unique to each mound, so just take min of the
%  interpolated surface.

for i = 1:numfiles
    gpsM = zeros(r(i).r.RasterSize);
    [xs, ys] = size(gpsM);
    gpsM = repmat(gpsM,1,1,3);                 
    stretchx = cos(ccent(i).cc(2)/3376200.00386);
    dx = [r(i).r.XWorldLimits(1):r(i).r.CellExtentInWorldX:r(i).r.XWorldLimits(2)];
    dy = [r(i).r.YWorldLimits(1):r(i).r.CellExtentInWorldY:r(i).r.YWorldLimits(2)];
    minz = min(min(a(i).a));
    for j = 1:length(dx)-1
        for k = 1:length(dy)-1
            gpsM(length(dy)-k,length(dx)-j,1) = real(sqrt((stretchx*(dx(j)-mat(i).mcent(1).mc(1)))^2+(dy(k)-mat(i).mcent(1).mc(2))^2));
            gpsM(length(dy)-k,length(dx)-j,2) = atan(sqrt((dy(k)-mat(i).mcent(1).mc(2))^2)/sqrt((dx(j)-mat(i).mcent(1).mc(1))^2));
            gpsM(k,j,3) = a(i).a(k,j)-minz;
        end
    end
    mat(i).mpolara(1).mp = gpsM;
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            for j = 1:length(dx)-1
                for k = 1:length(dy)-1
                    gpsM(length(dy)-k,length(dx)-j,1) = real(sqrt((stretchx*(dx(j)-mat(i).mcent(nmic).mc(1)))^2+(dy(k)-mat(i).mcent(nmic).mc(2))^2));
                    gpsM(length(dy)-k,length(dx)-j,2) = atan(sqrt((dy(k)-mat(i).mcent(nmic).mc(2))^2)/sqrt((dx(j)-mat(i).mcent(nmic).mc(1))^2));
                    gpsM(k,j,3) = a(i).a(k,j)-minz;
                end
            end
            mat(i).mpolara(nmic).mp = gpsM;
        end
    end
end

% 3376200.00386 meters per radian (or: horizontal distance * cos(YWorldLimits(1)/ <- ))