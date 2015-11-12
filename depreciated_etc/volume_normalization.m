vdistc = struct('vd',{[]});
ashaped = struct('as',{[]});
far = 0;

for i = 1:numfiles
    i
    mi = 1;
    ashaped(i).as = reshape(diffM(i).d,size(closest(i).c));
    mdsort = sortrows([closest(i).c', ashaped(i).as']);
    vdistc(i).vd(1) = 0;
    far = max(far, ceil(max(closest(i).c)));
        
    for rd = 2:ceil(max(closest(i).c))
        temp1 = 0;
        while mi < length(mdsort) && mdsort(mi,1) < rd
            temp1 = temp1 + max(mdsort(mi,2),0);
            mi = mi+1;
        end
        vdistc(i).vd(rd) = vdistc(i).vd(rd-1)+temp1;
    end
end

figure
hold on
for i = 1:numfiles
    plot(vdistc(i).vd ./ max(vdistc(i).vd))
end

figure
hold on
for i = 1:numfiles
    plot(vdistc(i).vd)
end