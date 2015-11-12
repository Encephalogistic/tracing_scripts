mdistc = struct('md',{[]});
cdistc = struct('cd',{[]});
%cudist = struct('c',{[]});
cudist2 = struct('c2',{[]});
metadist = [];
far = 0;

for i = 1:numfiles
    i
    mi = 1;
    ci = 1;
    mdsort = sort(closest(i).c);
    cdsort = sort(cclosest(i).cc);
    
    mdistc(i).md(1) = 0;
    cdistc(i).cd(1) = 0;
    cudist2(i).c2(1) = 0;
    
    if ceil(max(max(cclosest(i).cc,closest(i).c))) > far
        far = ceil(max(max(cclosest(i).cc,closest(i).c)));
    end
        
    for rd = 2:ceil(max(max(cclosest(i).cc,closest(i).c)))
        temp1 = 0;
        while mdsort(mi)<rd
            temp1 = temp1 + mdsort(mi);
            mi = mi + 1;
        end
        mdistc(i).md(rd) = mdistc(i).md(rd-1)+temp1;
                
        temp2 = 0;
        while cdsort(ci)<rd
            temp2 = temp2 + cdsort(ci);
            ci = ci + 1;
        end
        cdistc(i).cd(rd) = cdistc(i).cd(rd-1)+temp2;
               
        if temp2 == 0
            cudist2(i).c2(rd) = cudist2(i).c2(rd-1);
        else
            cudist2(i).c2(rd) = temp1/temp2 + cudist2(i).c2(rd-1);
        end
        %cudist2(i).c2 = cudist2(i).c2 ./ max(cudist2(i).c2);
    end
end

%    cudist(i).c = zeros([1,ceil(max(max(cclosest(i).cc,closest(i).c)))]);
%    for rd = ceil(r(i).r.CellExtentInWorldX/2):(ceil(max(cclosest(i).cc)) - (r(i).r.CellExtentInWorldX/2+1))
%        [mslope, yint] = polyfit([ceil(rd-r(i).r.CellExtentInWorldX/2):ceil(rd+r(i).r.CellExtentInWorldX/2)],...
%            mdistc(i).md(ceil(rd-r(i).r.CellExtentInWorldX/2):ceil(rd+r(i).r.CellExtentInWorldX/2)),1);
        
%        [cslope, cyint] = polyfit([ceil(rd-r(i).r.CellExtentInWorldX/2):ceil(rd+r(i).r.CellExtentInWorldX/2)],...
%            cdistc(i).cd(ceil(rd-r(i).r.CellExtentInWorldX/2):ceil(rd+r(i).r.CellExtentInWorldX/2)),1);
        
%        cudist(i).c(rd) = mslope(1) / cslope(1);

      %  cudist(i).c(rd) = mean(mdistc(i).md(ceil(rd-r(i).r.CellExtentInWorldX/2):ceil(rd+r(i).r.CellExtentInWorldX/2)))./...
        %    mean(cdistc(i).cd(ceil(rd-r(i).r.CellExtentInWorldX/2):ceil(rd+r(i).r.CellExtentInWorldX/2)));%+...
            %cudist(i).c(rd-1);
%    end
%    cudist(i).c(cudist(i).c == 0) = NaN;


figure
hold on
metadist = zeros(1,far);
for i = 1:numfiles
   cudist2(i).c2(isnan(cudist2(i).c2)) = max(cudist2(i).c2);
   cudist2(i).c2 = [cudist2(i).c2 ./ max(cudist2(i).c2), ones(1,length(metadist)-length(cudist2(i).c2))];
   plot(cudist2(i).c2)
   metadist = metadist + [cudist2(i).c2, zeros(1,length(metadist)-length(cudist2(i).c2))];
end
figure
plot(metadist./numfiles)