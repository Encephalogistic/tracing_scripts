%% Assign each point in the container (and mound) a distance to the nearest container edge.
%  Distance calculated using cross product with reference
%  to polygon vertices, adjusted for line segment termination points.  Also
%  plots naive cdf- but see mound_area_fraction for normalized graphs.

%% Major variable declarations

v1 = [];                                %List of vertices.
v2 = [];                                %Same list, offset by one.  Paired, defines a list of container edges.
closest = struct('c',{[]});             %Stores the shortest distance to any container wall for each mound point.
cclosest = struct('cc',{[]});           %Stores the shortest distance to any container wall for each container point.
gps = struct('g',{[]});                 %Links absolute matrix to raster values- gives GPS coordinates of each cell
                                        %(i.e. gps(i).g(5,10) gives X and Y
                                        % for a(i).a(5,10))
                 
%% Calculate the centroid of each container and mound.  Stolen without
%  remorse from HJ Sommer III (polygeom.m, see documentation).  Then,
%  assembles a vector assigning each point in the container and mound to a
%  distance value.  The vector doesn't record the original index in the
%  matrix, but preserves the left-to-right up-to-down reading order for
%  later reconstruction if desired.

for i = 1:numfiles
    ccent(i).cc = [];
    carea(i).ca = [];
    mcent(i).mc = [];
    marea(i).ma = [];
    
    contX = c(i).c.X(~isnan(c(i).c.X));
    contY = c(i).c.Y(~isnan(c(i).c.Y));
    mounX = m(i).m(1).X(~isnan(m(i).m(1).X));
    mounY = m(i).m(1).Y(~isnan(m(i).m(1).Y));
    
    [geom, centr, stuff] = polygeom(contX(1:end-1), contY(1:end-1));
    ccent(i).cc = [geom(2), geom(3)];
    carea(i).ca = geom(1);
    [geom, centr, stuff] = polygeom(mounX(1:end-1), mounY(1:end-1));
    mcent(i).mc = [geom(2) , geom(3)];
    marea(i).ma = geom(1);

    %Generate points to be measured.
    
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
    cclosest(i).cc = cclosest(i).cc .* reshape(double(a(i).a~=0),size(cclosest(i).cc));
    cclosest(i).cc(cclosest(i).cc == 0) = NaN;
    
    closest(i).c = cclosest(i).cc .* reshape(double(moundA(i).m~=0),size(cclosest(i).cc));
    closest(i).c(closest(i).c == 0) = NaN;

%    temp = cclosest(i).cc;

%    cclosest(i).cc(reshape(a(i).a~=0,size(cclosest(i).cc))) = 0;
%    cclosest(i).cc(reshape(logical(inContainer(i).ic),size(cclosest(i).cc)) &...
%    reshape(cclosest(i).cc == 0,size(cclosest(i).cc))) = temp;
    
%    closest(i).c = cclosest(i).cc .* reshape(double(moundA(i).m~=0),size(cclosest(i).cc));
 %   closest(i).c(reshape(logical(inContainer(i).ic),size(cclosest(i).cc)) && closest(i).cc == 0) = temp;
    
    figure
    hold on;
    surf(reshape(cclosest(i).cc,size(p)), 'EdgeColor','none');
end
%% Display naive cdf area plots

prodM = struct('p',{[]});
moundVector = [];
metaVector = [];
figure
hold on;

for i = 1:numfiles
    cdfplot(reshape(closest(i).c,1,[]))
    metaVector = [metaVector reshape(closest(i).c,1,[])];
end

figure
cdfplot(metaVector)