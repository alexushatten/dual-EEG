function Out = my_find_str(string,expr)

Out = find(~cellfun(@isempty,regexp(string,expr)));

end