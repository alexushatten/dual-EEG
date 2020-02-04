function DataStruct = robust_cell2struct(DataCell,Headers)
% convert cell array to structure using regular expressions to avoid problematic field names
%
% [DataStruct] = robust_cell2struct(DataCell, Headers)
%
%  Inputs
% --------
% DataCell: [n x n] cell array with table data (text or numbers)
% Headers: [n x 1 or 1 x n] cell array with headers of table data
%
%  Outputs
% ---------
% DataStruct: table data as structure with fields defined in Headers
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, February 2019
%-------------------------------------------------------------------------

DataStruct = struct;
for h = 1:length(Headers)
	DataStruct = setfield(DataStruct,regexprep(regexprep(regexprep(Headers{h},'\.|\/| |?|\(|\-|\)|:|''','_'),'é|è','e'),'à','a'),DataCell(:,h));
end

end
