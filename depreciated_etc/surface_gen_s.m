%%Reading the outputs from mound_cuts.m, interpolate the basal surface,
% sampling from the perimeter of the polygon and extrapolating inward.
% Should be formatted appropriately for later volume and relative distance
% calculations.  Switch commented lines for different point sampling.

%% Variable definitions.  
%  Also see parse_MOLA_mound_files.m and mound_cuts.m

index = struct('i',{[]});         % Stores selected a values for interpolating 
surfacecube = struct ('s',{[]});  % The interpolated surface, in a's coordinates

%% Sample points.
%  For current purposes, just populate the index with whatever onpolygon
%  happened to pick up.  Links those discrete coordinates to elevations.

for j = 1:numfiles
    
   % For 'random' onpoly sampling
   %lim = size(discX(j).onp);
   
   % For vertex-only sampling
   lim = size(discX(j).d(1:end-1));
   
   index(j).i = [];
   for k = 1:lim(2)
       
       % For 'random' onpoly sampling
       %index(j).i = [index(j).i a(j).a(discX(j).onp(k),discY(j).onp(k))];

       % For vertex-only sampling
       index(j).i = [index(j).i a(j).a( discX(j).d(k), discY(j).d(k))];
       
   end
end


%% Splining.
%  The cubic hermite spline is used here.  I think.  It looks nice and 
%  polynomial, anyway.

for j = 1:numfiles
    [xsize,ysize] = size(a(j).a);
    
    % For 'random' onpoly sampling
    %surfacecube(j).s = griddata(discX(j).onp,discY(j).onp,double(index(j).i),1:xsize,(1:ysize)', 'cubic'); 
    
    % For vertex-only sampling
    surfacecube(j).s = griddata(discX(j).d(1:end-1), discY(j).d(1:end-1),...
        double(index(j).i),1:xsize,(1:ysize)', 'cubic');
    
    surfacecube(j).s(isnan(surfacecube(j).s)) = 0; 
end