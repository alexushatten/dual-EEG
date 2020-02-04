
% Fitting parameters:

VerboseFiles = {'paths_to_VRB_verbose_files_from_best_grand_clustering_for_each_patient...'};

clear TemplateFile NumClusters
for v = 1:length(VerboseFiles)
    Vrb = read_vrb_file(VerboseFiles{v});
    TempLine = Vrb(~cellfun(@isempty,regexp(Vrb,'^Template file:')));
    if ispc
        TempLine2 = regexp(TempLine,':','split');
    else
        error('Cartool on UNIX... ha ha!')
    end
    TemplateFile{v} = [TempLine2{1}{2}(end),':',TempLine2{1}{3}];
    [~,FileName] = fileparts(TemplateFile{v});
    tempL = regexp(FileName,'\(..\)');
    tempR = regexp(FileName,'\(*\)');
    NumClusters(v) = str2double(FileName(tempL+1:tempR-1));
end
TemplateFile = TemplateFile';
NumClusters = NumClusters';

for v = 1:length(VerboseFiles)
    system(['explorer ',TemplateFile{v}]);
end

GCspikes_and_randomTemplates = {'path_to_EP_evoked_potential_files_of_templates_for_rest_from_best_grand_clustering_for_each_patient'};

CorrespondingSpikesTemplates = {'path_to_EP_evoked_potential_files_of_templates_for_spikes_from_best_grand_clustering_for_each_patient'};

TopoFiles = [GCspikes_and_randomTemplates;CorrespondingSpikesTemplates];

Topo = []; Nmaps = nan(length(TopoFiles),1); clear GroupNidx
for f = 1:length(TopoFiles)
    Temp = read_eph(TopoFiles{f});
    Topo = [Topo; Temp]; %#ok<AGROW>
    Nmaps(f) = size(Temp,1);
    GroupNidx{f} = (1:Nmaps(f))+sum(Nmaps(1:f-1)); %#ok<SAGROW>
end

Corr = corr(abs(Topo'));
Order = angular_order_eig(Corr);
figure; imagesc(Corr(Order,Order)); colorbar; title('Correlation between templates (sorted by angular order of eigenvectors)');

clear wPair VwPair
for f = 1:(length(TopoFiles)/2)
    [wPair{f},VwPair{f}] = map_r2c_mat(Corr,GroupNidx{f},GroupNidx{f});
end

wPair{1}






