v1 = [];
v2 = [];
gps = [0 0];
distM = struct('d',{[]});

%I have- diffM (with r), c.X, c.Y

for i = 1:numfiles
    v1 = [c(i).c.X(1:end-2) ; c(i).c.Y(1:end-2) ; zeros(1,length(c(i).c.Y(1:end-2)))];
    v2 = [v1(1:3,end), v1(1:3,1:end-1)];
    [xs,ys] = size(diffM(i).d);
    distM(i).d = zeros(size(diffM(i).d));
    
    for j = 1:xs
        for k = 1:ys
            [j, k, i]
            if (diffM(i).d(j,k) > 0)
                gps = [r(i).r.XWorldLimits(1)+j*r(i).r.CellExtentInWorldX r(i).r.YWorldLimits(1)+ ...
                    k*r(i).r.CellExtentInWorldY 0];
                for l = 1:length(v1)
                    a = v1(:,l)-v2(:,l);
                    b = gps' - v2(:,l);
                    distance = norm(cross(a,b))/norm(a);
                    if (distM(i).d(j,k) < distance)
                        distM(i).d(j,k) = distance;
                    end
                    %if (
                end
            end
        end
    end 
    end
    
% 3376200.00386 meters per radian (or: horizontal distance * cos(YWorldLimits(1)/ <- ))