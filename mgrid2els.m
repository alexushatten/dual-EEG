function mgrid2els(mgridFile)
%mgrid2els: writes .els (Cartool) based on .mgrid (BioImage Suite)
%
% .els will have the same basename and be written in folder containing
% .mgrid
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2018, updated October 2018
%-------------------------------------------------------------------------

ImDim = [256,256,256];

[elecMatrix, elecLabels]=mgrid2matlab_RM(mgridFile);

IsPartOfGrid = (~cellfun(@isempty,regexp(elecLabels,'[R|L]G_')))+1;
ContactNames = lower(regexprep(regexprep(regexprep(regexprep(regexprep(regexprep(regexprep(elecLabels,'LD_',''),'RD_',''),'RS_',''),'LS_',''),'RG_',''),'LG_',''),'_',''))';

% BioImage Suite is LIP whereas Cartool is RAS, and given Cartool XYZ and
% BioImage Suite IJK, if X <=> I, Y <=> K and Z <=> J
XYZ = [ImDim(1)-elecMatrix(:,1),ImDim(2)-elecMatrix(:,3),ImDim(3)-elecMatrix(:,2)];
    
write_els_file(spm_file(mgridFile,'ext','.els'),XYZ,ContactNames,0,IsPartOfGrid)



