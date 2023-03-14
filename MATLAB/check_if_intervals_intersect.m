function intersecting = check_if_intervals_intersect(lr1, lr2)
    diffs = lr1 - lr2'; % matrix, intersection happens if and only if there is a row or a column with a sign change
    if any( sign(diffs(:,1)) == -sign(diffs(:,2)))
        intersecting = true; return
    elseif any( sign(diffs(1,:)) == -sign(diffs(2,:)))
        intersecting = true;
    else
        intersecting = false;
    end
end