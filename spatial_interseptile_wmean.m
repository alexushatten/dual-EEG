function EEG = spatial_interseptile_wmean( EEG, D )
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
% EEG = spatial_interseptile_weighted_mean( EEG, D )
%
%  Inputs
% --------
% EEG: [channels x time] EEG traces
% D: distance matrix, e.g. result of:
%                         >> D = squareform(pdist(XYZ));
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

D(logical(eye(size(D))))=nan; % important for both Dinv and future search of neighbours
Dinv = 1./D; % inverse distance weight matrix
Dinv(logical(eye(size(D))))=1;

Neighbours = nan(size(D,1),Ix-1);
DistNeighbours = nan(size(D,1),Ix-1);
for c = 1:size(D,1)
    [SortedDist,SortedDistIdx] = sort(D(c,:));
    Neighbours(c,:) = SortedDistIdx(1:(Ix-1));
    DistNeighbours(c,:) = SortedDist(1:(Ix-1));
end

for t = 1:size(EEG,2)
    NewVals = zeros(size(EEG,1),1);
    for c = 1:size(EEG,1)
        Cluster = [c,Neighbours(c,:)];
        Seven = EEG(Cluster,t);
        MinSeven = Seven == min(Seven);
        MaxSeven = Seven == max(Seven);
        InInterseptile = (MinSeven + MaxSeven)<1;
        NewVals(c) = sum(vect(EEG(Cluster(InInterseptile),t))...
            .*vect(Dinv(Cluster(InInterseptile))))...
            /sum(Dinv(Cluster(InInterseptile)));
    end
    EEG(:,t) = NewVals;
end

end

