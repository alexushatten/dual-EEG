function Y = MAD_normalize_single_by_col( Y )
% => standardize matrix column by column  using the following formula : (X - mean(X)) / std(X)
% => if NaN values are ignored in the calculation of the mean and std but are
% not filtered
%
% Renaud Marquis @ LREN, 05.03.2015
% RM@FBMlab: update, May 2019, support single-precision float by converting to
% double...

if isa(Y,'single')
    Y = double(Y);
end
for i = 1:size(Y,2)
    Y(:,i) = (Y(:,i)-ignoreNaN(Y(:,i),@median))./ignoreNaN(Y(:,i),@mad); % normalize
end

end

