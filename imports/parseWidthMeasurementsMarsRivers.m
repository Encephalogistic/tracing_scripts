%parseWidthMeasurementsMarsRivers
%Edwin Kite
%Purpose: Take as input hand-picked bank pairs, and output paleo-channel
%width information.


%************************SET CONTROL PARAMETERS ***************************
if exist('WhichTransect','var') == 0;
   WhichTransect = 3; %Which geologic transect number
end

%Future improvement: replace fixed interp rms errors with dynamically loaded errors
%(using ArcGIS-output .TXT files).
if WhichTransect == 1
    interp_rms_error = [20.221 17.523 0 0 0];
elseif WhichTransect == 2
    interp_rms_error = [4.848 3.978 0 0 0]; %interpolation rms error
elseif WhichTransect == 3
    interp_rms_error = [67 67 0 0 0];
end
disp('Warning: Overriding interp rms with hand coded rms values from rms_data.txt files!')

distancecutoff = 500; %Distance, in m, above which line segments are no longer evaluated (for channel widths)
dottingstrategy = [1]; %if = 1, dot uniformly along channel banks; if = 0 (deprecated), find width from vertices.
        if sum(dottingstrategy==0)>=1 %Deprecated approach aggregating measurements at vertices only
            error('This dotting strategy is no longer implemented.')
        end
dbd = 2.5; %distance between dots (m)
Go8thMars2014Abstract = 0; %This includes binning routines as well as 8th-Mars-abstract specific figures
%************************SET GRAPHICS AND HOUSEKEEPING PARAMETERS *********
picturesque = 0;
picture_all_measurements = 0; %Time-intensive, used for debugging.
output_type = 'terse'; %'none', 'terse', or 'complete'
iss_for_display = 2; %Interpolated stratigraphic surface for display. 1 = planar; 2 = quadratic; 3 = Inverse-Distance Weighted; 4 = Kriged.
if iss_for_display == 4; error('Error bars for kriging not yet incorporated'); end

%************************IMPORT DATA FROM CSV FILES ***********************
loadWidthMeasurementsMarsRivers
%************************END OF IMPORT DATA FROM CSV FILES ****************

assert(min(px)==0); assert(min(py)==0); assert(minx~=0); assert(miny~=0);

%check image:
if picturesque >= 1
    figure
    scatter(px(zc(:,2)~=0),py(zc(:,2)~=0),10,z(zc(:,2)~=0)-zc(zc(:,2)~=0,2));
end    

%********* PEPPER POINTS UNIFORMLY ALONG CHANNEL BANKS **************
if sum(dottingstrategy==1)>=1
        ctl(1:max(bid)) = 0; %channel (actually bank) total length
        %Interpolate uniformly in s, dir coordinates
        for l = 1:(size(px,1)-1)
            if bid(l) == bid(l+1)%Is it from the same bank?
                lseg(l)     =   sqrt((px(l+1) - px(l)).^2 + ...
                                     (py(l+1) - py(l)).^2); %length of segment
                ctl(bid(l)) = ctl(bid(l)) + lseg(l); %bank total length
                lac(l+1)    = ctl(bid(l)); %length along bank for this vertex
            else %last point on channel trace
                %lac(l) = ctl(cid(l));
            end
        end
        
        assert(min(ctl(ctl>0))>dbd,'Error. Decrease distance between dots (dbd)')

        for m = 1:max(bid) 
            if sum(bid==m) >= 1
            %p = position vector of interpolated points
            p(m).x = interp1(lac(bid==m),px(bid==m),[0:dbd:ctl(m)]);
            p(m).y = interp1(lac(bid==m),py(bid==m),[0:dbd:ctl(m)]);
            p(m).z = interp1(lac(bid==m), z(bid==m),[0:dbd:ctl(m)]);
            if (WhichTransect == 1 || WhichTransect == 2 || WhichTransect == 3)
                for iss = 1:5 %iss = interpolated stratigraphic surface
                    p(m).zc(:,iss) = interp1(lac(bid==m), zc(bid==m,iss),[0:dbd:ctl(m)]);
                end
            end
            p(m).bid(1:max(size(p(m).x)))=m;
            p(m).bankunit(1:max(size(p(m).x)))=round(nanmean(bankunit(bid==m)));
            end
        end
end

%***************************FIND CHANNEL WIDTHS ***************************
dcls_temp = NaN;

for strategy_index = 1:length(dottingstrategy)
    dotuniformly = dottingstrategy(strategy_index)
    
    if dotuniformly == 0
        nstartpoints = size(px,1);
    elseif dotuniformly == 1
        clear udzc
        nstartpoints = max(size([p(:).bid]));
        udx = [p(:).x]; udy = [p(:).y]; udz = [p(:).z]; udbid = [p(:).bid]; %uniformly-dotted x, y, z, and bank-ID
         if (WhichTransect == 1 || WhichTransect == 2 || WhichTransect == 3)
            for iss = 1:5 %iss = interpolated stratigraphic surface
                 udzc = [p(1).zc];
                for b = 2:max([p(:).bid])
                 udzc = ...
                       [udzc;[p(b).zc]]; %concatenate vertically
                end
            end
                udzc = udzc';
         end
    end
    
%if dotuniformly == 0
for l = 1:nstartpoints %for each vertex
    if rem(l,100) == 0;disp(strcat('dotting strategy:',num2str(dotuniformly), ...
                                   ',l:',num2str(l)));end
    %if bid(min(size(px,1),l+1)) == bid(l) && bid(max(1,l-1)) == bid(l) ...
    %   && l > 1 && l < size(px,1)   %omit endpoints of banks
    
    if l > 1 && l < nstartpoints   %include endpoints of banks
    %find distance and coordinates of closest line segment
        if dotuniformly == 0
            x3 = px(l);y3 = py(l);
        elseif dotuniformly ==1
            x3 = udx(l);y3 = udy(l);
        end
        
        pd = sqrt((x3 - px).^2 + (y3 - py).^2);
        %pd(1) = distancecutoff*9e9; 
        pd(end) = distancecutoff*9e9; %prevent evaluation of first and last m-points
        xc = NaN*zeros(size(px));
        yc = xc;
        
        okmlist = find(pd<distancecutoff); %Do not consider very distant points as possible "closest" points.
        
    for mindex = 1:length(okmlist) 
           m = okmlist(mindex);
           %if pd(m) < distancecutoff
           %per http://stackoverflow.com/questions/6176227/for-a-point-in-an-irregular-polygon-what-is-the-most-efficient-way-to-select-th
           %but with sqrt added to denominator of equation for u(m) 
            x1 = px(m)  ; y1 = py(m);
            x2 = px(m+1); y2 = py(m+1);
        
        %segment is finite
            dQ1(m) = sqrt((x1 - x3)^2 + (y1 - y3)^2); %<-- Bugs result if these lines are moved inside the "filter out line segment between different bank ID's" loop
            
            dQ2(m) = sqrt((x2 - x3)^2 + (y2 - y3)^2);
            
        %DO consider the first line segment (for which bid(m-1)~=bid(m))
        if dotuniformly == 0
            acceptable_endpoint = (bid(m+1) == bid(m) && l~=m && bid(m)~=bid(l));
        elseif dotuniformly ==1
            acceptable_endpoint = (bid(m+1) == bid(m) && bid(m)~=udbid(l));
        end
            
        if acceptable_endpoint == 1 
      
%Following http://paulbourke.net/geometry/pointlineplane/
%using lambda instead of u and PC instead of P
lambda(m) = ( (x3-x1).*(x2-x1) + (y3-y1).*(y2-y1) ) ./ ...
             ((x2-x1).^2 + (y2-y1).^2);


if lambda(m)<0
    xc(m) = x1; yc(m) = y1;
    lambda(m) = 0;
elseif lambda(m)>1
    xc(m) = x2; yc(m) = y2;
    lambda(m) = 1;
else
    xc(m) = x1 + lambda(m)*(x2-x1);
    yc(m) = y1 + lambda(m)*(y2-y1);
end
          
        else
            %Skip. Zero values of xc will be used to set dcls_temp(m) to NaN
        end
            %Skip. Zero values of xc will be used to set dcls_temp(m) to NaN
    end
  
 dcls_temp = sqrt((x3 - xc).^2 + (y3 - yc).^2);        
 dcls_temp(isnan(xc)) = NaN; % Zero values of xc are used to set dcls_temp(m) to NaN
    
    [C,I] = nanmin(dcls_temp);
    if dotuniformly == 0
        lambdaforl(l) = lambda(I);
        clx(l) = xc(I);
        cly(l) = yc(I);
        dcls(l) = C;
        cls(l) = I; %line segment number, not point number
        clsbank(l) = bid(I);
    elseif dotuniformly == 1
        udlambdaforl(l) = lambda(I);
        udclx(l) = xc(I);
        udcly(l) = yc(I);
        uddcls(l) = C;
        udcls(l) = I; %line segment number, not point number
        udclsbank(l) = bid(I);
    end
    
     %We want to suppress the case where the closest *point* (not line segment) is a
            %vertex at the END of a bank (banks are indexed by 'bid')
                dQdebug(l) = 0;
            if I == 1; %The closest point is the first point in the list of all points on all banks
                dQdebug(l) = 9; %Suppress the width measurement
            elseif dQ1(I) == dcls_temp(I) && bid(I-1)~=bid(I)
                dQdebug(l) = 1;
            elseif dQ2(I) == dcls_temp(I) && bid(min(length(bid),I+2))~=bid(I)
                dQdebug(l) = 2;
            end
            
            if dQdebug(l) > 0
                if dotuniformly == 0
                    dcls(l) = NaN; cls(l) = NaN; clsbank(l) = NaN;
                elseif dotuniformly == 1
                    udcls(l) = NaN; uddcls(l) = NaN; udclsbank(l) = NaN;
                end
            end
    
    else
        if dotuniformly == 0
            cls(l) = NaN; dcls(l) = NaN; clsbank(l) = NaN;
        else
            udcls(l) = NaN; uddcls(l) = NaN; udclsbank(l) = NaN;
        end
    end
end
end

%******************* END OF FIND CHANNEL WIDTHS ***************************

if picturesque == 1
    figure %check plot
    scatter(px,py,'k*')
    hold on
    title('plot for debugging')
    EdwinColor = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};
    if picture_all_measurements == 1
            for l = 2:(size(px,1)-1)
                if isnan(cls(l)) == 0
                    line([px(l) clx(l)],[py(l) cly(l)],'Color', EdwinColor{1+rem(l,6)})
                    %line([px(l) px(cls(l))],[py(l) py(cls(l))],'Color', EdwinColor{1+rem(l,6)})
                    scatter(px(l),py(l),20,EdwinColor{1+rem(l,6)},'LineWidth',2)
                end
            end
            %    scatter(px(1:length(dQdebug)),py(1:length(dQdebug)),1500,dQdebug);

    for l = 2:max(size(udclx))
        if isnan(udcls(l)) == 0
            line([udx(l) udclx(l)],[udy(l) udcly(l)],'Color', EdwinColor{1+rem(l,6)})
            %line([px(l) px(cls(l))],[py(l) py(cls(l))],'Color', EdwinColor{1+rem(l,6)})
            scatter(udx(l),udy(l),20,EdwinColor{1+rem(l,6)},'LineWidth',2)
        end
    end

          %     scatter(px(1:length(dQdebug)),py(1:length(dQdebug)),50,dQdebug);
    end
end 
    
%*** AGGREGATE POINT-BY-POINT MEASUREMENTS TO FIND BANK-PAIR AVERAGES *****

%Find median width, dz, e.t.c. for each bank-pair
%The following comment is no longer relevant): Note that although we take log widths when aggregating the data and
%looking for trends with stratigraphic elevation, here we are only
%interested in (hopefully small) deviations from the true channel width.
%Therefore taking linear means is appropriate (and log vs. linear won't matter for the median in any case). 

%Each channel has two banks; they are assumed to be in order of this pairing!  
for nn = 1:(max(bid)/2) 
    css(nn)  = sum(bid==((nn*2) ) | bid==((nn*2)-1));%channel sample size (for calculating pooled variance, later)
    
    if css(nn)>0
        
        if sum(dottingstrategy==1)>=1 %Uniformly dotting points along channel banks
            s(11).cx(nn)         = nanmean(udx((udbid==(nn*2)        | udbid==((nn*2)-1)) & isnan(udcls)==0));
            s(11).cy(nn)         = nanmean(udy((udbid==(nn*2)        | udbid==((nn*2)-1)) & isnan(udcls)==0));
            s(11).cz(nn)         = nanmean(udz((udbid==(nn*2)        | udbid==((nn*2)-1)) & isnan(udcls)==0));
            s(11).czdev(nn)      =  nanstd(udz((udbid==(nn*2)        | udbid==(nn*2)-1) & isnan(udcls)==0));

            %pqual, pstyle and tienumber are the same for every vertex of a
            %bank-pair (including 'unpaired' vertices), so averaging them using bid or udbid returns the
            %same result. Therefore we use bid (bank ID) indexing and not
            %udbid (uniform-dotting bank ID) indexing.
            s(11).cqual(nn)      = nanmean(pqual((bid==(nn*2)     | bid==(nn*2)-1))); %note that FID == a(:,2)
            s(11).cpstyle(nn)    = nanmean(pstyle((bid==(nn*2)    | bid==(nn*2)-1)));
            if nanmean(tienumber((bid==(nn*2) | bid==(nn*2)-1))) ~= 0
                s(11).ctienumber(nn) = nanmean(tienumber((bid==(nn*2) | bid==(nn*2)-1)));
            else
                s(11).ctienumber(nn) = -nn; %Tied only to itself (negative tienumber)
            end
            
            if (WhichTransect == 1 || WhichTransect == 2 || WhichTransect == 3);
            for iss = 1:5 %channel, z (stratigraphic) elevation
                s(11).czs(nn,iss)   =  nanmean(udz((udbid==(nn*2) | udbid==(nn*2)-1)) ...
                    - udzc(iss,udbid==(nn*2) | udbid==(nn*2)-1)); 
                s(11).czsdev(nn,iss)=  nanstd( udz((udbid==(nn*2) | udbid==(nn*2)-1)) ...
                    - udzc(iss,udbid==(nn*2) | udbid==(nn*2)-1)); 
            end            
                s(11).czsdev(nn,:) = sqrt(s(11).czsdev(nn,:).^2 + interp_rms_error.^2);
            end 
            
            
            banktobankdistances = uddcls(udbid==(nn*2) | udbid==((nn*2)-1));
            s(11).cw(nn)      = nanmean(banktobankdistances);
            s(11).cwdev(nn)   = nanstd(uddcls(udbid==(nn*2) | udbid==(nn*2)-1));
            %Added 20 June 2014
            s(11).cwlog(nn)   = nanmean(log(banktobankdistances));
            s(11).cwdevlog(nn)= nanstd(log(banktobankdistances));
            
            s(11).cwrange(nn) = max(banktobankdistances)./min(banktobankdistances);
            s(11).cwmin(nn)   = nanmin(banktobankdistances);
            s(11).cwmax(nn)   = nanmax(banktobankdistances);
            %All points along a bank have the same 'bankunit' value
            s(11).cbankunit(nn) = nanmean(bankunit(bid==(nn*2) | bid==(nn*2)-1));
        end
                
        if exist('dz','var') == 1
            cdz(nn) = nanmean(dz((bid==(nn*2) | bid==(nn*2)-1) & isnan(dcls')==0));
        end
    else
        cz(nn) = NaN;
    end
end

%* END OF AGGREGATE POINT-BY-POINT MEASUREMENTS TO FIND BANK-PAIR AVERAGES*


%******** AGGREGATE BANK-PAIRS MEASURED ON THE SAME CHANNEL ***************

%Combine the mean and standard deviation of bank-pairs with the same
%tienumber to obtain a characteristic width for each tienumber.
%A common tienumber indicates that the bank-pairs probably sat on the same
%channel thread.

clear  bpti t
ti = 0;

for m = 1:length(s(11).cw) %For each channel
    if s(11).cqual(m)<8 %Do not include candidates
        
    if s(11).ctienumber(m) ~= -m %Negative tienumbers correspond to bank-pairs that do not have a tie  
        ctie = find(s(11).ctienumber == s(11).ctienumber(m));%bank-pairs with the same tienumber as this bank-pair
        ctie(s(11).cqual(ctie)>=8) = []; %If any of the tied bank-pairs are candidates, remove them from the
        %list of ties.
    else
        ctie = -m; %Tied only to itself (negative tienumber)
    end
    
    
    if (length(ctie)==1 || s(11).ctienumber(m) == -m) %This bank-pair has a unique tienumber.
        ti = ti + 1; %tie index
        t(ti).cw        = s(11).cw(m);
        t(ti).cwdev     = s(11).cwdev(m);
        %Added 20 June 2014.
        t(ti).cwlog     = s(11).cwlog(m);
        t(ti).cwdevlog  = s(11).cwdevlog(m);
        
        t(ti).cx        = s(11).cx(m);
        t(ti).cy        = s(11).cy(m);
        t(ti).cz        = s(11).cz(m);
        t(ti).czdev     = s(11).czdev(m);
        t(ti).cbankunit = s(11).cbankunit(m);
        bpti(m) = ti; %bank-pair tie index
        if (WhichTransect == 1 || WhichTransect == 2 || WhichTransect == 3); %stratigraphic elevations are available
           t(ti).czs    = s(11).czs(m,:);
           t(ti).czsdev = sqrt(s(11).czsdev(m,:).^2 + interp_rms_error.^2);
        end
        t(ti).intie = ctie;
        t(ti).tienumber = ctie;
    else %Combine the data for the bank-pairs with the same tienumber.
        if m == min(ctie); %only one t (tie) entry per tienumber
            
        ti = ti + 1; %tie index
        t(ti).cw        = nanmean(s(11).cw(ctie));
        %Added 20 June 2014
        t(ti).cwlog     = nanmean(s(11).cwlog(ctie));
        
        t(ti).cx        = nanmean(s(11).cx(ctie));
        t(ti).cy        = nanmean(s(11).cy(ctie));
        t(ti).cz        = nanmean(s(11).cz(ctie));
        t(ti).cbankunit = nanmean(s(11).cbankunit(ctie));
        t(ti).intie = ctie;    
        t(ti).tienumber = nanmean(s(11).ctienumber(ctie));
        if (WhichTransect == 1 || WhichTransect == 2 || WhichTransect == 3);
            t(ti).czs   = nanmean(s(11).czs(ctie,:));
        end
    %Based on Wikipedia article
    %'Standard_deviation#Combining_standard_deviations' as of 24FEB2014
    %Non-overlapping populations
    t(ti).cwdev = sqrt( ( sum(( [s(11).cwdev(ctie)].^2 + [s(11).cw(ctie)].^2 )) ...
                    / length(ctie))...
                    - t(ti).cw.^2);
    t(ti).czdev = sqrt( ( sum(( [s(11).czdev(ctie)].^2 + [s(11).cz(ctie)].^2 )) ...
                    / length(ctie))...
                    - t(ti).cz.^2);
    %Added 20 June 2014            
    t(ti).cwdevlog = sqrt( ( sum(( [s(11).cwdevlog(ctie)].^2 + [s(11).cwlog(ctie)].^2 )) ...
                    / length(ctie))...
                    - t(ti).cwlog.^2);            
                
    for iss = 1:5
    %Pool the variance for stratigraphic elevations for bank-pairs
    %with the same TieNumber. 'czsdev' here is the scatter in stratigraphic
    %elevations among the points dotted uniformly along each bank-pair.
    t(ti).czsdev(iss) = sqrt( ( sum(( [s(11).czsdev(ctie,iss)].^2 + [s(11).czs(ctie,iss)].^2 )) ...
                    / length(ctie))...
                    - t(ti).czs(iss).^2);
    %Now add RMS error for the interpolation as a whole. Except in the case of kriging error, 
    %this is uniform across the
    %interpolation surface.
    t(ti).czsdev(iss) = sqrt(t(ti).czsdev(iss).^2 + interp_rms_error(iss).^2);
    end              
        bpti(ctie) = ti;
        end
    end
    else
        bpti(m) = NaN;
    end
    
end

%Debug plots
if picturesque >= 1
figure
    errorbar([1:length(t)],[t.cw],[t.cwdev],'r.')
    hold on
    errorbar(0.25+bpti,[s(11).cw],[s(11).cwdev],'*');
    title('Widths, tied using tie-number vs individual bank-pairs')
figure    
    errorbar([1:length(t)],[t.cwlog],[t.cwdevlog],'r.')
    hold on
    errorbar(0.25+bpti,[s(11).cwlog],[s(11).cwdevlog],'*');
    title('log Widths, tied using tie-number vs individual bank-pairs')
        
figure
    errorbar(0.25+bpti,[s(11).cz],[s(11).czdev],'m*');
        hold on
    errorbar([1:length(t)],[t.cz],[t.czdev],'g.')
figure
    errorbar(0.25+bpti,[s(11).czs(:,iss_for_display)],[s(11).czsdev(:,iss_for_display)],'c*');
        hold on
    for ti = 1:length(t)
        errorbar([ti],[t(ti).czs(iss_for_display)],[t(ti).czsdev(iss_for_display)],'k.')
    end
end
    
%**************************** DISPLAY RESULTS *****************************

%Kludge fix for mixed preservation styles
s(11).cpstyle(s(11).cpstyle==9) = 8;
s(11).cpstyle(s(11).cpstyle==10) = 8;

if picturesque >= 1
figure%('units','normalized','outerposition',[0 0 1 1])
    h = errorbar(s(11).cz,s(11).cw,s(11).cwdev,'xk','LineWidth',2)
    hold on;grid on; colorbar; errorbar_tick(h,0)
    scatter(s(11).cz,s(11).cw,30,s(11).cx,'LineWidth',2)
        scatter(s(11).cz(s(11).cqual>=8),s(11).cw(s(11).cqual>=8),300,'ro','LineWidth',1)
    title('Raw elevation, candidates circled in red'); xlabel('Raw elevation (m)');ylabel('Channel width (m)');
    saveas(gca,strcat('figureWidthMeasurementsMarsRivers_Transect_',num2str(WhichTransect),'_1.eps'),'epsc');
    
%CORRECTED ELEVATION, INDIVIDUAL BANK-PAIRS:
figure%('units','normalized','outerposition',[0 0 1 1])
   %h = errorbar(cz(cqual<=4) - czc(cqual<=4,2)' ,s(11).cw(cqual<=4),s(11).cwdev(cqual<=4),'xk','LineWidth',1)
    h = errorbar(s(11).czs(s(11).cqual<=4,iss_for_display)' ,s(11).cw(s(11).cqual<=4),s(11).cwdev(s(11).cqual<=4),'xk','LineWidth',1)
    hold on;colorbar;grid on; errorbar_tick(h,0)
    set(gca,'FontSize',18)
    title('Stratigraphic elevation, color is log2(pstyle), no candidates, individual bank-pairs')
    xlabel('Stratigraphic elevation relative to F1-F2 contact (m)'); ylabel('Channel width (m)')
    %herrorbar(cz,cw,czdev,'LineSpec',' ') %scatter in NOMINAL ELEVATIONS
    %(not kridged stratigraphic position!) is usually negligible
    h = herrorbar(s(11).czs(s(11).cqual<=4,iss_for_display)',s(11).cw(s(11).cqual<=4),s(11).czsdev(s(11).cqual<=4,iss_for_display)')
    scatter(s(11).czs(s(11).cqual<=4,iss_for_display)' ,s(11).cw(s(11).cqual<=4),30,log(s(11).cpstyle(s(11).cqual<=4))./log(2),'LineWidth',2)    
    saveas(gca,strcat('figureWidthMeasurementsMarsRivers_Transect_',num2str(WhichTransect),'_2.eps'),'epsc');

    
%CORRECTED ELEVATION, ONE POINT PER TIE NUMBER:
   figure%('units','normalized','outerposition',[0 0 1 1])
        for irhm = 1:length(t)
            tzadjusted(irhm) = [t(irhm).czs(iss_for_display)];
            tzdevadjusted(irhm) = [t(irhm).czsdev(iss_for_display)];
        end
    h = errorbar(tzadjusted,[t(:).cw],[t(:).cwdev],'xk','LineWidth',1)
    hold on;grid on;colorbar; errorbar_tick(h,0);
    set(gca,'FontSize',18)
    title({'Stratigraphic elevation';'color is x position, one data point per TieNumber'})
    xlabel('Stratigraphic elevation relative to F1-F2 contact (m)'); ylabel('Channel width (m)')
    %herrorbar(cz,cw,czdev,'LineSpec',' ') %scatter in NOMINAL ELEVATIONS
    %(not kridged stratigraphic position!) is usually negligible
    h = herrorbar(tzadjusted,[t(:).cw],tzdevadjusted)
    scatter(tzadjusted,[t(:).cw],30,[t(:).cx],'LineWidth',2)
    saveas(gca,strcat('figureWidthMeasurementsMarsRivers_Transect_',num2str(WhichTransect),'_3.eps'),'epsc');
                                   
                    
end               
                    
%************************** SAVE OUTPUT ***********************************
switch output_type
    case 'terse'
        widthoutputname = strcat('outputWidthMeasurementsMarsRivers_Transect_',num2str(WhichTransect))        
        save(widthoutputname, ...
             's','t','minx','miny')
    case 'complete'
        disp('Not saving complete output.')
    case 'none'
end

%************************** END OF SAVE OUTPUT ****************************
                    
                                 
%******************************************For 8th Mars 2014 Abstract Submission *********************************                  
if Go8thMars2014Abstract == 1
figure('units','normalized','outerposition',[0 0 0.35 1])
    errorbar(s(11).cz(s(11).cqual<=4) - (cetl + cx(cqual<=4).*pfxc(WhichTransect) ...
                                      + cy(cqual<=4).*pfyc(WhichTransect)), ...
             s(11).cw(cqual<=4),s(11).cwdev(cqual<=4),'xr','LineWidth',2)
    hold on;grid on;
    czcorr = cz(cqual<=4) - (cetl + cx(cqual<=4).*pfxc(WhichTransect) ...
                                  + cy(cqual<=4).*pfyc(WhichTransect));
    set(gca,'FontSize',18)
    %title({'Planar correction, one point per pick on F1-F2 contact,';'color is log2(quality), no candidates'})
    xlabel('Stratigraphic elevation relative to F1-F2 contact (m)'); ylabel('Channel width (m)')
    %herrorbar(cz,cw,czdev,'LineSpec',' ') %scatter in NOMINAL ELEVATIONS
    %(not krided stratigraphic position!) is usually negligible
  %  scatter(cz(cqual<=4) - cetl + cx(cqual<=4).*pfxc(WhichTransect) ...
  %                               + cy(cqual<=4).*pfyc(WhichTransect)), ...
  %          s(11).cw(cqual<=4),30,log(cqual(cqual<=4))./log(2),'LineWidth',2)
    wunit0(1,1) =    nanmean(s(11).cw(czstrat <(0))); %Mean channel width above and below paleohydrologic break
    wunit0(2,1) =    nanmean(s(11).cw(czstrat>=(0))); %Mean channel width above and below paleohydrologic break
    wunit0(1,2) =  nanmedian(s(11).cw(czstrat <(0))); %Median channel width above and below paleohydrologic break
    wunit0(2,2) =  nanmedian(s(11).cw(czstrat>=(0))); %Median channel width above and below paleohydrologic break
    wunit0(1,3) =     nanstd(s(11).cw(czstrat <(0))); %s.d. of channel width above and below paleohydrologic break
    wunit0(2,3) =     nanstd(s(11).cw(czstrat>=(0))); %s.d. of channel width above and below paleohydrologic break
    
  line([0 0],[0.01 300],'Color','k')
  line([-150 0],[wunit0(1,1) wunit0(1,1)],'Color','k')
  line([0 150],[wunit0(2,1) wunit0(2,1)],'Color','k')
view([-90, 90])
set(gca,'ydir','reverse')
grid off;% set(gca,'yminortick','on')
set(gca,'yscale','log')
xlim([-70 130])
ylim([3 130]) 
set(gca,'ytick',[5 10 15 20 30 40 60 80 120])

                    saveas(gca,'For8thMars2014_DebugPlanarFit30MAR2014_4_PQualExcludeCandidates.fig','fig')
                    saveas(gca,'For8thMars2014_DebugPlanarFit30MAR2014_4_PQualExcludeCandidates.eps','epsc');

 %Split by x = 1e4 (for dtm seperation on transect 1:)
 
 figure('units','normalized','outerposition',[0 0 0.35 1])
    errorbar(cz(cqual<=4 & cx<1e4) - (cetl + cx(cqual<=4 & cx<1e4).*pfxc(WhichTransect) ...
                                  + cy(cqual<=4 & cx<1e4).*pfyc(WhichTransect)), ...
             s(11).cw(cqual<=4 & cx<1e4),s(11).cwdev(cqual<=4 & cx<1e4),'xr','LineWidth',2)
    hold on;grid on;
    czcorr = cz(cqual<=4 & cx<1e4) - (cetl + cx(cqual<=4& cx<1e4).*pfxc(WhichTransect) ...
                                  + cy(cqual<=4 & cx<1e4).*pfyc(WhichTransect));
    set(gca,'FontSize',18)
    title({'Planar correction, one point per pick on F1-F2 contact,';'color is log2(quality), no candidates'})
    xlabel('Stratigraphic elevation relative to F1-F2 contact (m)'); ylabel('Channel width (m)')
    %herrorbar(cz,cw,czdev,'LineSpec',' ') %scatter in NOMINAL ELEVATIONS
    %(not krided stratigraphic position!) is usually negligible
  %  scatter(cz(cqual<=4) - (cetl + cx(cqual<=4).*pfxc(WhichTransect) ...
  %                               + cy(cqual<=4).*pfyc(WhichTransect)), ...
  %          s(11).cw(cqual<=4),30,log(cqual(cqual<=4))./log(2),'LineWidth',2)
    wunit0(1,1) =    nanmean(s(11).cw(czstrat <(0) & cx<1e4)); %Mean channel width above and below paleohydrologic break
    wunit0(2,1) =    nanmean(s(11).cw(czstrat>=(0) & cx<1e4)); %Mean channel width above and below paleohydrologic break
    wunit0(1,2) =  nanmedian(s(11).cw(czstrat <(0) & cx<1e4)); %Median channel width above and below paleohydrologic break
    wunit0(2,2) =  nanmedian(s(11).cw(czstrat>=(0) & cx<1e4)); %Median channel width above and below paleohydrologic break
    wunit0(1,3) =     nanstd(s(11).cw(czstrat <(0) & cx<1e4)); %s.d. of channel width above and below paleohydrologic break
    wunit0(2,3) =     nanstd(s(11).cw(czstrat>=(0) & cx<1e4)); %s.d. of channel width above and below paleohydrologic break
    
  line([0 0],[0.01 300],'Color','k')
  line([-150 0],[wunit0(1,1) wunit0(1,1)],'Color','k')
  line([0 150],[wunit0(2,1) wunit0(2,1)],'Color','k')
view([-90, 90])
set(gca,'ydir','reverse')
grid off;% set(gca,'yminortick','on')
set(gca,'yscale','log')
xlim([-70 130])
ylim([3 130]) 
set(gca,'ytick',[5 10 15 20 30 40 60 80 120])

                    saveas(gca,'For8thMars2014_PQualExcludeCandidates_LowX.fig','fig')
                    saveas(gca,'For8thMars2014_PQualExcludeCandidates_LowX.eps','epsc');

%High split by x = 1e4

 figure('units','normalized','outerposition',[0 0 0.35 1])
    errorbar(cz(cqual<=4 & cx>1e4) - (cetl + cx(cqual<=4 & cx>1e4).*pfxc(WhichTransect) ...
                                  + cy(cqual<=4 & cx>1e4).*pfyc(WhichTransect)), ...
             s(11).cw(cqual<=4 & cx>1e4),s(11).cwdev(cqual<=4 & cx>1e4),'xr','LineWidth',2)
    hold on;grid on;
    czcorr = cz(cqual<=4 & cx>1e4) - (cetl + cx(cqual<=4& cx>1e4).*pfxc(WhichTransect) ...
                                  + cy(cqual<=4 & cx>1e4).*pfyc(WhichTransect));
    set(gca,'FontSize',18)
    title({'Planar correction, one point per pick on F1-F2 contact,';'color is log2(quality), no candidates'})
    xlabel('Stratigraphic elevation relative to F1-F2 contact (m)'); ylabel('Channel width (m)')
    %herrorbar(cz,cw,czdev,'LineSpec',' ') %scatter in NOMINAL ELEVATIONS
    %(not krided stratigraphic position!) is usually negligible
  %  scatter(cz(cqual<=4) - (cetl + cx(cqual<=4).*pfxc(WhichTransect) ...
  %                               + cy(cqual<=4).*pfyc(WhichTransect)), ...
  %          s(11).cw(cqual<=4),30,log(cqual(cqual<=4))./log(2),'LineWidth',2)
    wunit0(1,1) =    nanmean(s(11).cw(czstrat <(0) & cx>1e4)); %Mean channel width above and below paleohydrologic break
    wunit0(2,1) =    nanmean(s(11).cw(czstrat>=(0) & cx>1e4)); %Mean channel width above and below paleohydrologic break
    wunit0(1,2) =  nanmedian(s(11).cw(czstrat <(0) & cx>1e4)); %Median channel width above and below paleohydrologic break
    wunit0(2,2) =  nanmedian(s(11).cw(czstrat>=(0) & cx>1e4)); %Median channel width above and below paleohydrologic break
    wunit0(1,3) =     nanstd(s(11).cw(czstrat <(0) & cx>1e4)); %s.d. of channel width above and below paleohydrologic break
    wunit0(2,3) =     nanstd(s(11).cw(czstrat>=(0) & cx>1e4)); %s.d. of channel width above and below paleohydrologic break
    
  line([0 0],[0.01 300],'Color','k')
  line([-150 0],[wunit0(1,1) wunit0(1,1)],'Color','k')
  line([0 150],[wunit0(2,1) wunit0(2,1)],'Color','k')
view([-90, 90])
set(gca,'ydir','reverse')
grid off;% set(gca,'yminortick','on')
set(gca,'yscale','log')
xlim([-70 130])
ylim([3 130]) 
set(gca,'ytick',[5 10 15 20 30 40 60 80 120])

                    saveas(gca,'For8thMars2014_PQualExcludeCandidates_HighX.fig','fig')
                    saveas(gca,'For8thMars2014_PQualExcludeCandidates_HighX.eps','epsc');
                    
                    
                    
%END OF FOR 8TH MARS 2014                    

                    
figure
   
    for nn = 1:(max(bid)/2)
    pointwisezstrat = (zstrat((bid==(nn*2) | bid==(nn*2)-1)));
    plot(ones(size(pointwisezstrat))*(cz(nn) - (cetl + cx(nn).*pfxc(WhichTransect) ...
                                                     + cy(nn).*pfyc(WhichTransect))), ...
            pointwisezstrat);grid on; hold on
    end
        title('Check: should be a straight 1-to-1 line')
        ylabel('median elevation with stratigraphic correction');
        xlabel('pointwise median stratigraphic elevation');
   
        
figure
    title('Check')
    scatter(cx,cy,30,(cetl + cx.*pfxc(WhichTransect) ...
                           + cy.*pfyc(WhichTransect)),'*' ...
            );grid on;colorbar    
        hold on
    %scatter(cx,cy+100,30,-(czstrat - cz),'o' ...
    %        );grid on;colorbar    
    %    hold on       
    xlabel('x');ylabel('y')
figure

    errorbar(cdz,s(11).cw,cwdev,'xk','LineWidth',2)
    hold on
    scatter(cdz,s(11).cw,50,'co','fill')
    scatter(cdz(cdz<0),         s(11).cw(cdz<0),50,'bo','fill')
    scatter(cdz(cdz>0 & cdz<50),s(11).cw(cdz>0 & cdz<50),50,'wo','fill')
    scatter(cdz(cdz>0 & cdz<50),s(11).cw(cdz>0 & cdz<50),50,'ko')

    %errorbar(cdz(cwhichdtm==2),cw(cwhichdtm==2),cwdev(cwhichdtm==2),'xr')

    %herrorbar(cdz,cw,crootsemi,'x')
xlabel('Relative Time (Height Above F1/F2 Contact (m))','FontSize',20)
ylabel('River Channel Width (m)','FontSize',20)
set(gca,'FontSize',20)
set(gca,'xminortick','on')
set(gca,'yminortick','on')
ylim([0 70])
xlim([-100 150])
%grid on
box on
line([-100,0],[29.3,29.3],'Color','b','LineWidth',2)
line([50,150],[15.8616,15.8616],'Color','c','LineWidth',2)
line([0,50],[ 26.8001, 26.8001],'Color','k','LineWidth',2)


%error('Not doing probability summation while debugging')

%Quick summation of probability distribution (assuming uncorrelated, Guassian error
%bars) From Wikipedia, "Multivariate normal distribution"
p = zeros(601,401);
warning('Assuming fixed sigmax for debug')
for nn = 1:(max(bid)/2) %add the probability kernel of each observation
oldp = p;
nn
            sigmax = 29.0847;%<-- Assuming fixed sigma x for debugging %crootsemi(nn);
            %The value of fixed sigma x comes from
            %f1_f2_contact_second_pass_b20g02_only_plane_fit_one_point_per_contact_pick_24feb14_rms_output
            
            sigmay = s(11).cwdev(nn);
            normfactor = 1/(2*pi*sigmax*sigmay);
            mux = cz(nn) - (cetl + cx(nn).*pfxc(WhichTransect) ...
                                 + cy(nn).*pfyc(WhichTransect));
            muy = s(11).cw(nn);
            
            if sum(isnan([sigmax,sigmay,mux,muy]))==0
                for x = -300:1:300 %height above F1/F2 contact (m)
                    for y = -100:1:300 %river channel width (m)

                            p(x+301,y+101) = p(x+301,y+101) + ...
                            normfactor*exp(-1/2*(((x-mux)^2/(sigmax^2))+((y-muy)^2/(sigmay^2))));
                    end
                end
            end
                %For debug only
                %check that the additional probability sums to 1 (note that
                %the discretisation and the cutoffs will introduce errors!)
                %if abs((sum(p(:)) - sum(oldp(:)))-1) > 0.01
                %    error('MVN Dist Checksum Error!')
                %end
 end

 [X,Y] = meshgrid(-300:1:300,-100:1:300);
%contour(X,Y,p');grid on %show the sum of all kernels

%show 1 sigma and 2 sigma error windows for channel width as a function of
%stratigraphic position

pcs = cumsum(p,2);
ps = sum(p');
for l = 1:601 %cumsum normalized
pcsn(l,:) = pcs(l,:)./ps(l);
end

hold on
contour(X,Y,pcsn',[0.05,(0.31731/2),1-((0.31731/2)),0.95]);

%find pooled variance for 50m stratigaphic-position bins: (Note to self:
%should preferably use *fractional* variance?)
stratbin = 1;
for stratz = -200:50:200
    to_use(stratbin).v = find(cdz>stratz & cdz <stratz+50);
    ninbin(stratbin) = length(to_use(stratbin).v);
    pooledvar(stratbin) = sqrt(sum(css(to_use(stratbin).v).*(cwdev(to_use(stratbin).v).^2))/sum(css(to_use(stratbin).v)));
    stratbin = stratbin + 1;
end

%Overall pooled variance:
opooledvar = sqrt(sum(css.*(cwdev.^2))/sum(css))

end
