function [LessPairable, ValLessPairable, Mapping, Values] = find_less_pairable_mat_element( Mat, RowIdx, ColIdx )
% find_less_pairable_mat_element: match rows to columns of matrix such that
% pairs are ranked using highest value (e.g. highest correlation) and identify
% elements - if any - in rows or columns which cannot be mapped
%
% [LessPairable, ValLessPairable, Mapping, Values] = ...
%					find_less_pairable_mat_element( Mat, RowIdx, ColIdx )
%
%  Inputs
% --------
% Mat: [n x m] symmetric matrix (e.g. correlation matrix)
% RowIdx: [j x 1 or 1 x j] indices of rows to consider in Mat
% ColIdx: [k x 1 or 1 x k] indices of rows to consider in Mat
%         (There should be no overlap between ColIdx and RowIdx except in
%         the case where ColIdx is strictly equivalent to RowIdx, in which
%         case the matrix is assumed to be symmetric)
%
%  Outputs
% ---------
% LessPairable: less pairable pair of elements in Mat within RowIdx and ColIdx
% ValLessPairable: value in Mat corresponding to LessPairable
% Mapping: [p x 2] array of correspondances between RowIdx and ColIdx,
%          where p = min([j,k])
% Values: [p x 1] vector of values of Mat associated to each pair in
%         Mapping. When the matrix is symmetric and has an odd number of
%         rows (and columns), the last row that cannot be matched to any
%         column (except to itself) has a NaN value.
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018, NOT FULLY TESTED /!\
%-------------------------------------------------------------------------

[Mapping, Values] = map_r2c_mat( Mat, RowIdx, ColIdx );
MatNaN = Mat;
MatNaN(logical(eye(size(MatNaN))))=nan;

if any(RowIdx(Mapping(:,1))==ColIdx(Mapping(:,2))) % checks if symmetric
    RemainingElements = find(RowIdx(Mapping(:,1))==ColIdx(Mapping(:,2))); % if such then the last element is the only one remaining
    ValLessPairable = nan; %#ok<NASGU>
end
if length(RowIdx)~=length(ColIdx)
    if length(RowIdx)>length(ColIdx) % more rows than columns
        RemainingElements = setdiff(RowIdx,RowIdx(Mapping(:,1)));
    elseif length(RowIdx)<length(ColIdx) % more columns than rows
        RemainingElements = setdiff(ColIdx,ColIdx(Mapping(:,2)));
    end
    if length(RemainingElements)>1 % if more than 1 remaining element
        warning('EXPERIMENTAL --- be careful with results interpretation')
        if length(RowIdx)>length(ColIdx) % nrows>ncols
            Mapping2 = map_r2c_mat_nonsym( MatNaN(RemainingElements,ColIdx), 1:length(RemainingElements), 1:length(ColIdx));
%             LessPairable = [RemainingElements(Mapping2(:,1))',Mapping2(:,1)];
            LessPairable = [RemainingElements(Mapping2(:,1))',nan(size(Mapping2,1),1)];
        elseif length(RowIdx)<length(ColIdx) % ncols>nrows
            Mapping2 = map_r2c_mat_nonsym( MatNaN(RowIdx,RemainingElements), 1:length(RowIdx), 1:length(RemainingElements));
%             LessPairable = [Mapping2(:,1),RemainingElements(Mapping2(:,2))'];
            LessPairable = [nan(size(Mapping2,1),1),RemainingElements(Mapping2(:,2))'];
        end
        ValLessPairable = nan;
        if size(Mapping2,1)<length(RemainingElements)
%             warning('Your are passing a special case to the function (matrix size has one dimension which is more than the double of the other and remaining elements do not disappear after one iteration), which will use an experimental subfunction to try to handle your situation... Be careful with the results !')
%             MatBKP = Mat; Mat = Mat';
%             [RowIdx, ColIdx] = swap(RowIdx, ColIdx);
%             [LessPairable, ~, Mapping] = find_less_pairable_mat_element( Mat, ColIdx, RowIdx);
%             LessPairable = fliplr(LessPairable);
%             ValLessPairable = MatBKP(LessPairable(1),LessPairable(2));
%             Mapping = fliplr(Mapping);
%             Values = MatBKP(Mapping);
            error('Case is not supported')
%             if length(RowIdx)>length(ColIdx) % nrows>ncols
%                 RemainingElements2 = setdiff(RowIdx,RowIdx(Mapping2(:,1)));
%                 [LessPairable, ValLessPairable, Mapping, Values] = flpme_subfun( MatNaN(RemainingElements2,ColIdx), 1:length(RemainingElements2), 1:length(ColIdx) );
%             elseif length(RowIdx)<length(ColIdx) % ncols>nrows
%                 RemainingElements2 = setdiff(ColIdx,ColIdx(Mapping2(:,2)));
%                 [LessPairable, ValLessPairable, Mapping, Values] = flpme_subfun( MatNaN(RowIdx,RemainingElements2), 1:length(RowIdx), 1:length(RemainingElements2) );
%             end
        end
    else
        if length(RowIdx)>length(ColIdx) % more rows than columns
            RemainingElements = find(RemainingElements==RowIdx);
        elseif length(RowIdx)<length(ColIdx) % more columns than rows
            RemainingElements = find(RemainingElements==ColIdx);
        end
        if length(RowIdx)>length(ColIdx)
%             [~,J] = max(MatNaN(RemainingElements,:));
            LessPairable = [RemainingElements,nan];
        elseif length(RowIdx)<length(ColIdx)
%             [~,J] = max(MatNaN(:,RemainingElements));
            LessPairable = [nan,RemainingElements];
        end
        ValLessPairable = nan;%MatNaN(LessPairable(1),LessPairable(2));
    end
else
    LessPairable = Mapping(end,:);
    ValLessPairable = Values(end);
end
    
end

% function [LessPairable, ValLessPairable, Mapping, Values] = flpme_subfun( Mat, RowIdx, ColIdx )
% 
% warning('Your are passing a special case to the function (matrix size has one dimension which is more than the double of the other and remaining elements do not disappear after one iteration), which will use an experimental subfunction to try to handle your situation... Be careful with the results !')
% 
% [Mapping, Values] = map_r2c_mat_nonsym( Mat, RowIdx, ColIdx );
% MatNaN = Mat;
% MatNaN(logical(eye(size(MatNaN))))=nan;
% 
% if any(RowIdx(Mapping(:,1))==ColIdx(Mapping(:,2))) % checks if symmetric
%     RemainingElements = find(RowIdx(Mapping(:,1))==ColIdx(Mapping(:,2))); % if such then the last element is the only one remaining
%     ValLessPairable = nan;
% end
% if length(RowIdx)~=length(ColIdx)
%     if length(RowIdx)>length(ColIdx) % more rows than columns
%         RemainingElements = setdiff(RowIdx,RowIdx(Mapping(:,1)));
%     elseif length(RowIdx)<length(ColIdx) % more columns than rows
%         RemainingElements = setdiff(ColIdx,ColIdx(Mapping(:,2)));
%     end
%     if length(RemainingElements)>1 % if more than 1 remaining element
%         if length(RowIdx)>length(ColIdx) % nrows>ncols
%             MatNaN(RemainingElements,ColIdx)
%             [Mapping2, ValLessPairable] = map_r2c_mat_nonsym( MatNaN(RemainingElements,ColIdx), 1:length(RemainingElements), 1:length(ColIdx));
%             LessPairable = [RemainingElements(Mapping2(:,1))',Mapping2(:,1)];
%         elseif length(RowIdx)<length(ColIdx) % ncols>nrows
%             MatNaN(RowIdx,RemainingElements)
%             [Mapping2, ValLessPairable] = map_r2c_mat_nonsym( MatNaN(RowIdx,RemainingElements), 1:length(RowIdx), 1:length(RemainingElements));
%             LessPairable = [Mapping2(:,1),RemainingElements(Mapping2(:,2))'];
%         end
%         if size(Mapping2,1)<length(RemainingElements)
%             if length(RowIdx)>length(ColIdx) % nrows>ncols
%                 [LessPairable, ValLessPairable, Mapping, Values] = flpme_subfun( MatNaN(RemainingElements,ColIdx), 1:length(RemainingElements), 1:length(ColIdx) );
%             elseif length(RowIdx)<length(ColIdx) % ncols>nrows
%                 [LessPairable, ValLessPairable, Mapping, Values] = flpme_subfun( MatNaN(RowIdx,RemainingElements), 1:length(RowIdx), 1:length(RemainingElements) );
%             end
%         end
%     else
%         if length(RowIdx)>length(ColIdx)
%             [~,J] = max(MatNaN(RemainingElements,:));
%             LessPairable = [RemainingElements,J];
%         elseif length(RowIdx)<length(ColIdx)
%             [~,J] = max(MatNaN(:,RemainingElements));
%             LessPairable = [J,RemainingElements];
%         end
%         ValLessPairable = MatNaN(LessPairable(1),LessPairable(2));
%     end
% else
%     LessPairable = Mapping(end,:);
%     ValLessPairable = Values(end);
% end
% 
% end