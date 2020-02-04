function X = empty_cells_to_nan_array(X)

EmptyX = cellfun(@isempty,X);

X(EmptyX) = {nan};

try
    X = cell2mat(X);
catch
    X = cellfun(@double,X);
end


end