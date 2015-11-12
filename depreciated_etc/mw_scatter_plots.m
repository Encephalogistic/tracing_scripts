moundWidth = struct('mw',{[]});
containerArea = struct('ca',{[]});
containerRelief = struct('cr',{[]});
moundRelief = struct('mr',{[]});

for i = 1:numfiles
    moundW = 0;
    index = 1;
    
    while moundW == 0
        if cudist2(i).c2(index) >= 0.5
            moundW = index
        end
        index = index + 1;
    end
    moundWidth(i).mw = moundW;
end

for i = 1:numfiles
    containerArea(i).ca = sum(sum(~isnan(a(i).a./a(i).a)))
end

for i = 1:numfiles
    tempInterp = wideInterp(i).w;
    tempInterp(tempInterp  == 0) = NaN;
    containerRelief(i).cr = sqrt((max(max(tempInterp))-min(min(tempInterp)))^2)
end

for i = 1:numfiles
    moundRelief(i).mr = max(max(diffM(i).d))
end

figure
hold on
for i = 1:numfiles
    if i == 5 ||...
        i == 6 ||...
        i == 7 ||...
        i == 9||...
        i == 10||...
        i == 12||...
        i == 14
        scatter(containerArea(i).ca,moundWidth(i).mw,100,'red','filled')
    elseif i == 4||...
            i == 23
        scatter(containerArea(i).ca,moundWidth(i).mw,100,'green','filled')
    else
        scatter(containerArea(i).ca,moundWidth(i).mw,100,'blue','filled')
    end
    title = ('Moat width as a function of container area')
end

figure
hold on
for i = 1:numfiles
    if i == 5 ||...
        i == 6 ||...
        i == 7 ||...
        i == 9||...
        i == 10||...
        i == 12||...
        i == 14
        scatter(containerRelief(i).cr,moundWidth(i).mw,100,'red','filled')
    elseif i == 4||...
            i == 23
        scatter(containerRelief(i).cr,moundWidth(i).mw,100,'green','filled')
    else
        scatter(containerRelief(i).cr,moundWidth(i).mw,100,'blue','filled')
    end
    title = ('Moat width as a function of container relief')
end

figure
hold on
for i = 1:numfiles
    if i == 5 ||...
        i == 6 ||...
        i == 7 ||...
        i == 9||...
        i == 10||...
        i == 12||...
        i == 14
        scatter(moundRelief(i).mr,moundWidth(i).mw,100,'red','filled')
    elseif i == 4||...
            i == 23
        scatter(moundRelief(i).mr,moundWidth(i).mw,100,'green','filled')
    else
        scatter(moundRelief(i).mr,moundWidth(i).mw,100,'blue','filled')
    end
    title = ('Moat width as a function of mound relief')
end