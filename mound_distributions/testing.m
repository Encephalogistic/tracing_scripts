M = [1 2 3;1 2 3];
fid = fopen('output.csv','w');
fprintf(fid,'%s\r\n',['One,Two,Three']);
fclose(fid);
dlmwrite('output.csv', M,'-append','delimiter',',','roffset',0,'coffset',0);

N = csvread('output.csv',1,0)

    match = false;
    if ~isempty(allFiles)
        for j = 1:length(allFiles(:,1))
            if (i == allFiles(j,1) && shapename == allFiles(j,2))
                match = true
            end
        end
    end
    if match == false
        allFiles = [allFiles;xmeters(i),ymeters(i),dist(i),angle(i),distc(i),thick(i),elev(i),slope(i),azimuth(i),...
            sllb(i),slub(i),azlb(i),azub(i)];
    end
    
    if isempty(allFiles(allFiles(:,1) == i-1 && allFiles(:,2) == shapename))
    end