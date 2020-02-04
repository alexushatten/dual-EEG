function ProximCheck = dual_get_close_mrk( Check, EOF, Lag )
% dual_get_close_mrk: get markers close in time
%
% ProximCheck = dual_get_close_mrk( Check, EOF, Lag )
%
%  Inputs
% --------
% Check: cell array, output of dual_check_mrk.m
% EOF : integer, end of signals (normally both T and Y should end at same time!)
% Lag : allowed lag (in samples) between timings of T and Y.
%
%  Outputs
% ---------
% ProximCheck: structure, with pairs of markers, number of times they are
% close according to Lag
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, August 2018
%-------------------------------------------------------------------------

warning off %#ok<WNOFF>
Comb = nchoosek(1:size(Check,1),2);
% Comb = [Comb,fliplr(Comb)];
% => #RM commented line above because we only care about the "Hits" here so we don't need it in both directions
Hits = zeros(size(Comb,1),1);
for n = 1:size(Comb,1)
    prc_for_loop(n,size(Comb,1),1);
    for f = 1:numel(Check{1,3}) % for each .mrk file
        HitsTemp = timings_matching(...
        [Check{Comb(n,1),3}{f},Check{Comb(n,1),4}{f}],...
        [Check{Comb(n,2),3}{f},Check{Comb(n,2),4}{f}],...
        EOF,Lag);
        if isnan(HitsTemp)
            HitsTemp = 0;
        end
        Hits(n) = Hits(n)+HitsTemp;
    end
end
ProximCheck = [Check(Comb(find(Hits),1),1),Check(Comb(find(Hits),2),1),cellstr(num2str(Hits(find(Hits))))]; %#ok<FNDSB> % #RM: we need find because "Hits" can be greater than 1

warning on %#ok<WNON>

end
