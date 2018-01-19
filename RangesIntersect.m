function intersect_vals = RangesIntersect(range1, range2,range3)
% range1 = errorz(1,:);
% range2 = errorz(2,:);
% range3 = errorz(3,:);
if range1(1) < range2(1)
    lowRange1 = range1;
    highRange1 = range2;
else
    lowRange1 = range2;
    highRange1 = range1;
end
if range2(1) < range3(1)
    lowRange2 = range2;
    highRange2 = range3;
else
    lowRange2 = range3;
    highRange2 = range2;
end
if range1(1) < range3(1)
    lowRange3 = range1;
    highRange3 = range3;
else
    lowRange3 = range3;
    highRange3 = range1;
end
intersect_vals = zeros(1,3);
intersect_vals(1) = lowRange1(2)>highRange1(1);
intersect_vals(2) = lowRange2(2)>highRange2(1);
intersect_vals(3) = lowRange3(2)>highRange3(1);