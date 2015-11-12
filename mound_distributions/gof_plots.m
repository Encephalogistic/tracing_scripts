wid = [];

for i = 1:numfiles
    wid = [wid mat(i).moundWidth(1).mw * bin];
    for nmic = 2:10
        if length(m(i).m)>=nmic
            wid = [wid mat(i).moundWidth(nmic).mw * bin];
        end
    end
end

bins = 1:8:max(wid);
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

[logBinCounts] = histc(log(wid),1:.5:10);
[logNData,logNGOF] = fit((1:length(logBinCounts))',logBinCounts','gauss1');
figure
plot(logNData,1:length(logBinCounts),logBinCounts)

%% Chi2 goodness of fit.

[yn pval] = chi2gof(wid,'nbins',9);
normal_chi2_yn = yn
normal_chi2_pval = pval

[yn pval] = chi2gof(log(wid),'nbins',9);
lognormal_chi2_yn = yn
normal_chi2_pval = pval

figure
hwid = histogram(wid,9)
hstats = hwid.BinWidth*.5 + hwid.BinEdges(1:end-1);
n_obs = hwid.Values;

n_exp = normpdf(hstats,mean(wid),std(wid));
n_exp = n_exp/sum(n_exp);
n_exp = sum(n_obs)*n_exp;

subplot(1,2,1),bar(hstats,n_obs,'r')
subplot(1,2,2),bar(hstats,n_exp,'b')

chi2calc_normal = sum((n_obs - n_exp).^2 ./ n_exp)
chi2crit_normal = chi2inv(.95,6)

figure
hwid = histogram(log(wid),9)
hstats = hwid.BinWidth*.5 + hwid.BinEdges(1:end-1);
n_obs = hwid.Values;

n_exp = normpdf(hstats,mean(log(wid)),std(log(wid)));
n_exp = n_exp/sum(n_exp);
n_exp = sum(n_obs)*n_exp;

subplot(1,2,1),bar(hstats,n_obs,'r')
subplot(1,2,2),bar(hstats,n_exp,'b')

chi2calc_lognorm = sum((n_obs - n_exp).^2 ./ n_exp)
chi2crit_lognorm = chi2inv(.95,6)

figure
hwid = histogram(wid,6)
hstats = hwid.BinWidth*.5 + hwid.BinEdges(1:end-1);
n_obs = hwid.Values;

n_exp = exppdf(hstats,mean(wid));
n_exp = n_exp/sum(n_exp);
n_exp = sum(n_obs)*n_exp;

subplot(1,2,1),bar(hstats,n_obs,'r')
subplot(1,2,2),bar(hstats,n_exp,'b')

chi2calc_exp = sum((n_obs - n_exp).^2 ./ n_exp)
chi2crit_exp = chi2inv(.95,7)

