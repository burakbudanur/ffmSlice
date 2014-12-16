function sc = SliceCond(x)
%SliceCond(x) = 0 if x is on the slice hyperplane

T = Lg();
xp = slicep();
tp = T*xp;

sc = tp'*x;
