%parse_HiRISE tiff files and do downslope things

numf = 15; 
ftan(1).f = [];
meangrad(1).m = [];
notm(1).n = [];
shslope(1).ss = 0;
shazimuth(1).sa = 0;
elslope(1).es = 0;
elazimuth(1).ea = 0;
vecelev(1).ve = [];
container = 8;

for l = 1:numf
    % Read in elevation .tif files
    l
    name1745 = ['D:\jsneed_work\gale_sulfate_canyon\ESP_012907_1745_ESP_013540_1745\isolated_trace_tiffs\ESP_012907_1745_ESP_013540_1745_align_1-DEM-adj_' int2str(l-1) '.shp.tif'];
    name1745
    [elev(l).e,ras(l).rs] = geotiffread(name1745);
    elev(l).e(elev(l).e<-2e4) = NaN;
    terrelev = elev(l).e;
    
    % Create a best-fit bedding plane for the surface elevation map
    [p,q] = meshgrid(linspace(ras(l).rs.XWorldLimits(1),ras(l).rs.XWorldLimits(2),ras(l).rs.RasterSize(2)),...
        linspace(ras(l).rs.YWorldLimits(2),ras(l).rs.YWorldLimits(1),ras(l).rs.RasterSize(1)));
    [elslope(l).es,elazimuth(l).az,er,ste,norm(l).n] = regressplanemod(p(:),q(:),elev(l).e(:));
       
    % Find the steepest angle for elevation map and bedding surface
    [fx, fy] = gradient(elev(l).e);
    fx = mean(mean(fx(~isnan(fx))))
    fy = mean(mean(fy(~isnan(fy))))
    
    [gx, gy] = gradient(BVec(l).B(1)*p + BVec(l).B(2)*q + BVec(l).B(3));
    gx = mean(mean(gx(~isnan(gx))))
    gy = mean(mean(gy(~isnan(gy))))
    
    % Generate characteristic points along best-fit planes, and find
    % relevant angles between them
    X1(l) = p(ceil(length(p(:,1))/2),ceil(length(p(1,:))/2));
    X2(l) = p(ceil(length(p(:,1))/2),ceil(length(p(1,:))/2  + fx*100));
    X3(l) = p(ceil(length(p(:,1))/2),ceil(length(p(1,:))/2  + gx*100));
    Y1(l) = q(ceil(length(q(:,1))/2),ceil(length(q(1,:))/2));
    Y2(l) = q(ceil(length(q(:,1))/2 + fy*100),ceil(length(q(1,:))/2));
    Y3(l) = q(ceil(length(q(:,1))/2 + gy*100),ceil(length(q(1,:))/2));
    
    RP1(l) = BVec(l).B(1)*X1(l) + BVec(l).B(2)*Y1(l) + BVec(l).B(3);
    RP2(l) = BVec(l).B(1)*X2(l) + BVec(l).B(2)*Y2(l) + BVec(l).B(3);
    RP3(l) = BVec(l).B(1)*X3(l) + BVec(l).B(2)*Y3(l) + BVec(l).B(3);
    RQ1(l) = norm(l).n(1)*X1(l) + norm(l).n(2)*Y1(l) + norm(l).n(3);
    RQ2(l) = norm(l).n(1)*X2(l) + norm(l).n(2)*Y2(l) + norm(l).n(3);
    
    PA(l) = rad2deg(asin((RP1(l)-RP2(l))/sqrt((X1(l)-X2(l))^2 + (Y1(l)-Y2(l))^2 + (RP1(l)-RP2(l))^2)));
    QA(l) = rad2deg(asin((RQ1(l)-RQ2(l))/sqrt((X1(l)-X2(l))^2 + (Y1(l)-Y2(l))^2 + (RQ1(l)-RQ2(l))^2)));
    
    angdif(l) = PA(l)-QA(l)
    
    %Plot 3d map of surface, sampled points, best-fit plane of bedding
    %surface, and angle vectors (blue = elevation azimuth, green = bedding
    %plane azimuth, red = downslope component of bedding surface)
    figure
    hold on
    surf(ras(l).rs.XLimWorld(1)+[1:ras(l).rs.XLimIntrinsic(2)]*ras(l).rs.DeltaX',ras(l).rs.YLimWorld(2)+ ...
        [1:ras(l).rs.YLimIntrinsic(2)]*ras(l).rs.DeltaY',elev(l).e,'EdgeColor','none')
    quiver3(X1(l),Y1(l),RQ1(l),X1(l)-X2(l),Y1(l)-Y2(l),RQ1(l)-RQ2(l),'b')
    surf(ras(l).rs.XLimWorld(1)+[1:ras(l).rs.XLimIntrinsic(2)]*ras(l).rs.DeltaX',ras(l).rs.YLimWorld(2)+ ...
        [1:ras(l).rs.YLimIntrinsic(2)]*ras(l).rs.DeltaY',BVec(l).B(1)*p + BVec(l).B(2)*q + BVec(l).B(3),'EdgeColor','none')
    quiver3(X1(l),Y1(l),RQ1(l),X1(l)-X2(l),Y1(l)-Y2(l),RP1(l)-RP2(l),'r')
    quiver3(X1(l),Y1(l),RQ1(l),X1(l)-X3(l),Y1(l)-Y3(l),RP1(l)-RP3(l),'g')
    scatter3(divS(l).mx,divS(l).my,divS(l).ras,'g','filled')
    title(['Trace ',num2str(l)])
    
end

%Plot gradient vectors of each trace, spatially.  (Full HiRISE image is not
%imported in to Matlab, so you just have to imagine a nice contour map.)
figure
hold on
for i = 1:numf
    quiver(X1(i),Y1(i),X1(i)-X2(i),Y1(i)-Y2(i),10,'b')
    quiver(X1(i),Y1(i),X1(i)-X3(i),Y1(i)-Y3(i),10,'r')
    text(X1(i) + 50,Y1(i) + 50,int2str(i))
end