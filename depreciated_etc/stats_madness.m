
wid = [];

for i = 1:numfiles
    wid = [wid mat(i).moundWidth(1).mw];
    for nmic = 2:10
        if length(m(i).m)>=nmic
            wid = [wid mat(i).moundWidth(nmic).mw];
        end
    end
end

bins = 1:8:165;
[binCounts] = histc(wid,bins);

[normalData,normalGOF] = fit((1:length(binCounts))',binCounts','gauss1');
figure
plot(normalData,1:length(binCounts),binCounts)

[powerData,powerGOF] = fit((1:length(binCounts))',binCounts','power1');
figure
plot(powerData,1:length(binCounts),binCounts)

[expData,expGOF] = fit((1:length(binCounts))',binCounts','exp1');
figure
plot(expData,1:length(binCounts),binCounts)

[logBinCounts] = histc(log(wid),1:.25:10);
[logNData,logNGOF] = fit((1:length(logBinCounts))',logBinCounts','gauss1');
figure
plot(logNData,1:length(logBinCounts),logBinCounts)