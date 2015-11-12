prodM = struct('p',{[]});
moundVector = [];
metaVector = [];
hold on;

for i = 1:numfiles
    distM(i).d(distM(i).d == 0) = NaN;
    prodM(i).p = distM(i).d;
    moundVector = reshape(prodM(i).p,1,[]);
    cdfplot(moundVector);
    metaVector = [metaVector moundVector];
end

cdfplot(metaVector)