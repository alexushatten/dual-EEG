function EEG = spatial_interseptile_weighted_mean( EEG, XYZ )
% spatial_interseptile_weighted_mean: spatial filter as performed by
% Cartool software, instantaneous filter which removes local outliers by
% spatially smoothing the maps without losing its topographical
% characteristics: for each electrode, the values of the 6 closest
% neighbours are determined, plus the central electrode value itself; the 7
% data points are sorted; the minimal and maximal values are removed by
% dropping the first and last items of this list; the remaining values are
% then averaged, with weights proportional to the inverse distance to the
% central electrode. The central electrode is given a weight of 1.
%
% EEG = spatial_interseptile_weighted_mean( EEG, XYZ )
%
%  Inputs
% --------
% EEG: [channels x time] EEG traces
% XYZ: [channels x 3] array of coordinates (these should NOT be scaled /
%           normalized / etc. ! (like e.g. for spline interpolation)
%
%  Outputs
% ---------
% EEG: spatially filtered EEG
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

Ix = 7; % you can change the number of neighbours considered (and the type of inter-xxx that will be used, e.g. interquartile, interseptile, interdecile, ...)

D = squareform(pdist(XYZ)); % euclidean distance

if ispc % #RM@FBMlab: #TODO: compile MEX for UNIX as well
    if size(EEG,1)==size(D,1) && size(D,1) == size(D,2)
        EEG = spatial_interseptile_wmean_mex(EEG,D);
    end
else
    % #RM@FBMlab: might be interesting to try out polar distance later...
    D(logical(eye(size(D))))=nan; % important for both Dinv and future search of neighbours
    Dinv = 1./D; % inverse distance weight matrix
    Dinv(logical(eye(size(D))))=1;
    % Dinv = scale_mat_0_1(Dinv); % does not solve the problem, because then
    % the weight for the central electrode is too low...

    Neighbours = nan(size(XYZ,1),Ix-1);
    DistNeighbours = nan(size(XYZ,1),Ix-1);
    for c = 1:size(XYZ,1)
        [SortedDist,SortedDistIdx] = sort(D(c,:));
        Neighbours(c,:) = SortedDistIdx(1:(Ix-1));
        DistNeighbours(c,:) = SortedDist(1:(Ix-1));
    end

    for t = 1:size(EEG,2)
        prc_for_loop(t,size(EEG,2),100);
        NewVals = zeros(size(EEG,1),1);
        for c = 1:size(EEG,1)
            Cluster = [c,Neighbours(c,:)];
            Seven = EEG(Cluster,t);
            MinSeven = Seven == min(Seven);
            MaxSeven = Seven == max(Seven);
    %         OutInterseptile = (MinSeven + MaxSeven)>0;
            InInterseptile = (MinSeven + MaxSeven)<1;
            NewVals(c) = sum(vect(EEG(Cluster(InInterseptile),t))...
                .*vect(Dinv(Cluster(InInterseptile))))...
                /sum(Dinv(Cluster(InInterseptile)));
        end
        EEG(:,t) = NewVals;
    end
end

end

