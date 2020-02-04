
% Some electrodes should not be corrected because the electrode has
% been taken out of the stylus and does not follow a straight line
% despite being a depth electrode... !

Paths2mgrid = {
    'E:\DATA\sub-50.mgrid'
    };

for sub = 1:length(Paths2mgrid)
    
    mgridFname = Paths2mgrid{sub};
    
    [eCoords, elecLabels, elecRgb, elecPairs, elecPresent]=mgrid2matlab_RM(mgridFname);
    % eCoords=eCoords-1; % Make coordinates same as in mgrid file (thus first slice has a coordinate of 0, last has a coordinate of 255)
    % (as in makeIniLocTxtFile.m)
    
    % Get groups of electrodes
    [Bipoles,ElecMatch,labelsB,labelsB2] = bipolar_montage(lower(regexprep(regexprep(elecLabels,'[L|R]D',''),'_','')'));
            
    % ASSUMPTION: each electrode has its own unique label after underscore ("xxx" in pattern [R|L][G|S|D]_xxx)
    % do it per electrode, not for all at a time!!
    NcontactsElec = [0;cumsum(cellfun(@length,match_vectors(unique(ElecMatch,'stable'),ElecMatch)))]+1;
    
%     if sub == 8
%         % find(~cellfun(@isempty,regexp(elecLabels,'LD_IAG*')))
%         % Fortunately these are the 8 first contacts, it will easier to
%         % fix...
%         NcontactsElec = NcontactsElec(2:end);
%         AlignEqCoords = eCoords(1:8,:);
%         
%     else
    
        AlignEqCoords = [];
    
%     end
    
    for n = 2:(length(NcontactsElec))
        
        % ASSUMPTION: cordinates have to be sorted along each dimension (n => 1) !
        AlignEqCoords = [AlignEqCoords; AlignEqualizeContactDepthElec(eCoords(NcontactsElec(n-1):NcontactsElec(n)-1,:))]; %#ok<AGROW>
        
    end
    
%     copyfile(mgridFname,spm_file(mgridFname,'suffix','_align_equ'))
    % Write new .mgrid file
    matlab2mgrid_RM(mgridFname,AlignEqCoords,1);
    
end

% Write also .els file for visualization in Cartool:
for sub = 1:length(Paths2mgrid)
    
    mgridFname = Paths2mgrid{sub};
    mgrid2els(spm_file(mgridFname,'suffix','_align_equ'));

end
