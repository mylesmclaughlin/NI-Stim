%--------------------------------------------------------------------------
function PhaseDiff_removed = removejumps(PhaseDiff, tol)
if nargin <2
    tol = pi;
end
PhaseDiff_removed = PhaseDiff;
for r = 1:length(PhaseDiff_removed)
    if PhaseDiff_removed(r)< -pi
        PhaseDiff_removed(r)= PhaseDiff_removed(r)+2*pi;
    elseif PhaseDiff_removed(r)> pi
        PhaseDiff_removed(r) = PhaseDiff_removed(r)-2*pi;
    end
end