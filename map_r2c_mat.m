function [Mapping, Values] = map_r2c_mat( Mat, RowIdx, ColIdx )
% map_r2c_mat: maps rows to columns based on highest value in array in a 1
% to 1 fashion (i.e. avoids a row to be mapped to multiple columns)
%
% [ Mapping, Values ] = map_r2c_mat( Mat, RowIdx, ColIdx )
%
%  Inputs
% --------
% Mat: [n x m] array (e.g. correlation matrix)
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

SymFLAG = false;
if isempty(setdiff(RowIdx,ColIdx))
    [nn,mm] = size(Mat);
    if nn~=mm
        error('Matrix is not square although indices are the same, exiting...')
    else
        SymFLAG = true;
    end
    fprintf('Indices are the same, assuming symmetric matrix...\n')
    if ~issymmetric(round(Mat,10))
        warning('Matrix does not seem to be symmetric at 10 decimals precision!')
    end
    if all(diag(Mat)==1)
        %         Mat(logical(~triu(ones((nn)),1)))=nan;
        Mat(logical(eye(nn)))=nan;
    else
        warning('Diagonal of square matrix is not equal to 1, assuming it has been already NaNed, make sure it is the case!')
    end
else
    if ~all(setdiff(RowIdx,ColIdx)==RowIdx) || ~all(setdiff(ColIdx,RowIdx)==ColIdx)
        error('Indices are (only) partially overlapping, check input arguments')
    end
end

MatRanked = sort(vectorize(Mat(RowIdx,ColIdx)),'descend');
MatRanked = MatRanked(~isnan(MatRanked));
if numel(MatRanked)~=numel(unique(MatRanked))
    % if such, first fix original matrix, then re-create MatRanked,
    % otherwise there might be no match anymore later...
    % ... warn user:
    warning('Elements of input matrix are not unique although it is not symmetric, fix for this situation is to substract small constant to duplicate element in an arbitray way !');
    if mean(RowIdx)>mean(ColIdx)
        Mask = tril(ones(size(Mat)),1);
    else
        Mask = triu(ones(size(Mat)),1);
    end
    MatTriangleUp = Mat(logical(Mask));
    while numel(MatTriangleUp)~=numel(unique(MatTriangleUp))
        TempDup = find(cellfun(@length,match_vectors(MatTriangleUp,MatTriangleUp,1))>1);
        warning('Maximum error induced by correction: %d',length(TempDup)*eps)
        MatTriangleUp(TempDup) = MatTriangleUp(TempDup) + ((1:length(TempDup)).*eps)';
    end
    
    Mat(logical(Mask)) = MatTriangleUp;
    Mat(~logical(Mask)) = nan;
    % re-create:
    MatRanked = sort(vectorize(Mat(RowIdx,ColIdx)),'descend');
    MatRanked = MatRanked(~isnan(MatRanked));
end

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
        if ~SymFLAG
            [ii(v),jj(v)] = find(Mat(RowIdx,ColIdx)==MatRanked(w)); %#ok<AGROW>
        else
            % if is is odd and symmetric matrix, one of the element cannot be
            % paired:
            if v==length(RowIdx) && logical(rem(length(RowIdx),2))
                TheLastOne = setdiff(1:length(RowIdx),ii(1:v-1));
                % % or equivalently:
                % TheLastOne = setdiff(ColIdx,jj(v));
                ii(v) = TheLastOne; %#ok<AGROW>
                jj(v) = TheLastOne; %#ok<AGROW>
            else
                [ii(v),jj(v)] = find(Mat(RowIdx,ColIdx)==MatRanked(w),1); %#ok<AGROW>
            end
        end
        if ~SymFLAG
            if length(find(ii(v)==ii(1:v)))<2 && length(find(jj(v)==jj(1:v)))<2
                v = v+1; w = w+1;
            else
                w= w+1;
            end
        else
            if ~(v==length(RowIdx) && logical(rem(length(RowIdx),2)))
                if length(find(ii(v)==ii(1:v)))<2 && length(find(jj(v)==jj(1:v)))<2
                    ii(v+1) = jj(v); %#ok<AGROW>
                    jj(v+1) = ii(v); %#ok<AGROW>
                    v = v+2; w = w+2;
                else
                    w= w+1;
                end
            end
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
