%%Takes the generated .mat file, and extracts the points that occupy mound
% outlines (discrete coordinates, for matrix manipulation).  Prep for
% origin assignment, basal surface interpolation, and volume calculations.
% MOLA_mounds_working.mat must be preloaded

%%TODO: expand to handle muliple mounds per container

%% Variables, defined.  Because I usually code Java, that's why.
%  For .mat variables, see parse_MOLA_mound_files.m

discX = struct('d',{[]},'onp',{[]});    % ordered x-coordinates for the polygon vertices and edges
discY = struct('d',{[]},'onp',{[]});    % ordered y-coordinates for the polygon vertices and edges
moundA = struct('m',{[]});              % new matrix generated zeroing the exterior of a mound shape
inContainer = struct('ic',{[]});       % boolean matrix specifying whether a point is in container i

% smaller matrix limited to polygon dimensions (calculation efficiency
% measure)
maxX = 0;
maxY = 0;
minX = 0;
minY = 0;
moundAred = struct('m',{[]}, 'maxx', maxX, 'maxy', maxY, 'minx', minX, 'miny', minY);



%% Recast shapefile info in discrete coordinates, and find 
% polygon maximums/minimums appropriately

for i = 1:numfiles
    discX(i).d = zeros(size(m(i).m(1).X));
    discY(i).d = zeros(size(m(i).m(1).Y));
   [discX(i).d,discY(i).d] = worldToDiscrete(r(i).r,m(i).m(1).X,m(i).m(1).Y); 
   moundAred(i).maxx = max(discX(i).d);
   moundAred(i).maxy = max(discY(i).d);
   moundAred(i).minx = min(discX(i).d);
   moundAred(i).miny = min(discY(i).d);
end

%% Generate the moundA and moundAred matrices as zeros, then populate it with values only
% from within the given mound shape, leaving everything else as 0.
% Populates the ordered list of points on the polygon vertices for later,
% nefarious uses.  (i.e. surface_gen files)

for i = 1:numfiles
   i
   moundA(i).m = zeros(size(a(i).a));
   inContainer(i).ic = zeros(size(a(i).a));
   moundAred(i).m = zeros(moundAred(i).maxx+2-moundAred(i).minx,moundAred(i).maxy+2-moundAred(i).miny);
   [e,f]=size(moundAred(i).m);
   
   for j = 1:e
       for k = 1:f
           [inp, onp] = inpolygon(j-1+moundAred(i).minx,k-1+moundAred(i).miny,discX(i).d,discY(i).d);
           if inp || onp
             moundAred(i).m(j,k) = a(i).a(j-1+moundAred(i).minx,k-1+moundAred(i).miny);
             moundA(i).m(j-1+moundAred(i).minx,k-1+moundAred(i).miny) =...
                 a(i).a(j-1+moundAred(i).minx,k-1+moundAred(i).miny);
             inContainer(i).ic(j-1+moundAred(i).minx,k-1+moundAred(i).miny) = 1;
           end
           
           if onp
             discX(i).onp = [discX(i).onp j-1+moundAred(i).minx];
             discY(i).onp = [discY(i).onp k-1+moundAred(i).miny];
           end
       end
   end
end