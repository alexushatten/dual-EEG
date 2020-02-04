function Timings = prune_timings( Timings, TimeDiff )
% prune_timings: prune timings that are too close in time based on a
% time difference threshold using recursive removal
%
% Timings = prune_timings( Timings, TimeDiff )
%
%  Inputs
% --------
% Timings: [n x 1] numeric vector of timings / indices
%
% TimeDiff: [1 x 1] numeric, time difference threshold
%
%  Outputs
% ---------
% Timings: pruned Timings, separated by at least TimeDiff
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

% Timings with time difference below threshold are given by
% find(diff(Timings)>TimeDiff), but because there might be a serie of
% timings that are too close in time, e.g. 3rd, 4th, 5th, 6th, as soon as
% we remove one of these, the next one will be fine. We don't want to
% remove all timings that follow the first one, but rather spare only the
% ones that are sufficiently separated, so we will recursively remove only
% the first one: 
if any(diff(Timings)<TimeDiff)
    RemoveMe = find(diff(Timings)<TimeDiff,1)+1;
    Timings(RemoveMe) = [];
end

% Recursively call this function until no timing is below TimeDiff
if any(diff(Timings)<TimeDiff)
    Timings = prune_timings( Timings, TimeDiff );
end

end

