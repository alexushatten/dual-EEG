function [Mapping, Values] = map_r2c_mat_nonsym( Mat, RowIdx, ColIdx )
% map_r2c_mat: maps rows to columns based on highest value in array in a 1
% to 1 fashion (i.e. avoids a row to be mapped to multiple columns)
%
% /!\ USE map_r2c_mat.m RATHER THAN THIS FUNCTION, the latter is only for
% special cases and should be called by find_less_pairable_mat_elements.m !
%
% [ Mapping, Values ] = map_r2c_mat_nonsym( Mat, RowIdx, ColIdx )
%
%  Inputs
% --------
% Mat: [n x m] non-symmetric array (e.g. not a correlation matrix!)
% RowIdx: [j x 1 or 1 x j] indices of rows to consider in Mat
% ColIdx: [k x 1 or 1 x k] indices of rows to consider in Mat
%         (There should be no overlap between ColIdx and RowIdx except in
%         the case where ColIdx is strictly equivalent to RowIdx, in which
%         case the matrix is assumed to be symmetric)
%
%  Outputs
% ---------
% Mapping: [p x 2] array of correspondances between RowIdx and ColIdx,
%          where p = min([j,k])
% Values: [p x 1] vector of values of Mat associated to each pair in
%         Mapping. When the matrix is symmetric and has an odd number of
%         rows (and columns), the last row that cannot be matched to any
%         column (except to itself) has a NaN value.
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

MatRanked = sort(vectorize(Mat(RowIdx,ColIdx)),'descend');
MatRanked = MatRanked(~isnan(MatRanked));

w=1;v=1; % iter=0;
Temp = {};

[~,MaxDim] = max([length(RowIdx),length(ColIdx)]);
if MaxDim == 1 %&& MaxDim > length(RowIdx)
    V4temp = 1:length(RowIdx);
else
    V4temp = 1:length(ColIdx);
end

% "chaise musicale darwinienne"
while iscell(Temp) && w < length(MatRanked) && v <= length(RowIdx) && v <= length(ColIdx) % iter<10000
    %     iter = iter+1;
        
    [ii(v),jj(v)] = find(Mat(RowIdx,ColIdx)==MatRanked(w)); %#ok<AGROW>
    
    
    if length(find(ii(v)==ii(1:v)))<2 && length(find(jj(v)==jj(1:v)))<2
        v = v+1; w = w+1;
    else
        w= w+1;
    end
    
    if MaxDim==1
        Temp = match_vectors(V4temp,ii,1);
    else
        Temp = match_vectors(V4temp,jj,1);
    end
end

Mapping = [ii',jj'];

Values = zeros(size(Mapping,1),1);
for k = 1:length(ii)
    Values(k) = Mat(RowIdx(ii(k)),ColIdx(jj(k)));
end

end
