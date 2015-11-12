%% Integrates over a mound volume from the upper bound of observed
%  surface height to the lower bound of interpolated basal surface.
%  IBS may be imported as either a cubic spline or a kriging.
%  Also preserves upper-lower at each discrete cell in the matrix.

%% Define Variables
%  Also see parse_MOLA_mound_files.m, mound_cuts.m, surface_gen_s, etc.

wideInterp = struct('w',{[]}); 
diffM = struct('d',{[]});
moundVolume = struct('m',0);

%% Create a donut, then fill it.

for i = 1:numfiles
    inverseA = double(moundA(i).m == 0) .* double(a(i).a);
    wideInterp(i).w = inverseA;
    inZone = inverseA == 0;
    wideInterp(i).w = wideInterp(i).w + inZone .* surfacecube(i).s';
end

%% Create a new matrix of the difference between observed subaerial surface
%  and interpolated basal surface at each matrix coordinate.

for i = 1:numfiles
    diffM(i).d = double(a(i).a) - double(wideInterp(i).w);
    diffM(i).d(diffM(i).d == 0) = NaN;
    moundVolume(i).m = sum(sum(diffM(i).d));
end
%% Produce the volume of each mound.  Trapezoidal approximation is preferred.
%  Horizontal resolution = 300m