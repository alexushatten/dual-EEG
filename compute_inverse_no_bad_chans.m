function [ RISvect, RISnorm, RISsvd, OptimalReg ] = compute_inverse_no_bad_chans( Data, ISfile, Regularization )
% compute_inverse_no_bad_chans: compute results of inverse solutions given
% EEG traces and inverse solution matrix (.is file from Cartool), ignoring
% channels marked as bad.
%
% [ RISvect, RISnorm, RISsvd, OptimalReg ] = compute_inverse_no_bad_chans( Data, ISfile )
% [ RISvect, RISnorm, RISsvd, OptimalReg ] = compute_inverse_no_bad_chans( Data, ISfile, Regularization )
%
%  Inputs
% --------
% Data: structure with fields:
%                   - EEG: [channels x time] EEG traces
%                   - bad: [n x 1] vector with indices of bad channels
% ISfile: char, path to inverse solution file (output of Cartool)
% Regularization (optional): integer, Tikhonov regularization to use
%       (remember that regularization "0" should be requested with "1",
%       and instead of 0:12 we have 1:13 in Matlab!)...
%       Alternatively, Regularization can be set to 'optimal' (default),
%       which will calculate the optimal regularization based on L-corner
%       of the Tikhonov regularization values against the average norm of
%       all solution points.
%
%  Outputs
% ---------
% RISvect: results of inverse solution (product of EEG and inverse
%          solution), in vectorial form (intensity of dipole in x, y and z
%          directions)
% RISnorm: same as RISvect but taking the norm of each dipole, hence
%          reducing dimensionality
% RISsvd: same as RISvect but taking the first eigenvector (SVD) of each
%         dipole, hence reducing dimensionality, as with VOI extraction as
%         done in SPM and described in details in Rubega et al (2019)
%         (doi.org/10.1007/s10548-018-0691-2)
% OptimalReg: index of optimal regularization (based on L-corner)
%
% See also COMPUTE_INVERSE
%
%-------------------------------------------------------------------------
% NB: to get SVD across each ROI, use RISvect and compute SVD externally,
%     as here it will be done on the base of ROIs definition
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2019, last edited November 2019
%-------------------------------------------------------------------------

% Read inverse solution:
[IShdr,IS] = fbmlab_read_is(ISfile);
ISr = reshape(IS,3,IShdr.numsolutionpoints,IShdr.numelectrodes,IShdr.numregularizations);
clear IS % free memory

if size(ISr,3)~=size(Data.EEG,1)
    ISr(:,:,Data.bad,:) = [];
end

if nargin < 3
    Output = 'all';
    warning('Requesting sources for all regularization values can produce huge arrays and result in out-of-memory errors! You have been warned...')
    if nargout > 3
        error('Please input Regularization as ''optimal'' if you want the optimal regularization index')
    end
else
    if isnumeric(Regularization)
        Output = 'one';
        if nargout > 3
            error('Cannot output optimal regularization if it is given as input')
        end
    elseif strcmpi(Regularization,'optimal')
        Output = Regularization;
        if nargout < 4
            warning('Optimal regularization requested but not third output argument, results of inverse solutions for all regularization values will be outputed but optimal one will not be accessible')
        end
    else
        error('Invalid third input argument %s, please specify index or ''optimal''',Regularization)
    end
end

FlagHuge = false;
if strcmpi(Output,'all') || strcmpi(Output,'optimal')
    if strcmpi(Output,'all')
        RISvect = cell(IShdr.numregularizations,1);
    end
    fprintf('Computing sources in vectorial form...\n')
    AvgNorm = nan(IShdr.numregularizations,1);

	% A simple product will be calculated. The matrices might be huge, however.
    % Even if the matrices can hold in RAM, the speed of subsequent calculations
    % can be greatly accelerated if the objects are converted to single-precision float.
    % Unless some values are very small and might lead to loss of precision, we convert
    % to single-precision by default:
    AllMax = [max(abs(ISr(:))),max(abs(Data.EEG(:)))];
    AllMin = [min(abs(ISr(:))),min(abs(Data.EEG(:)))];
    if ~(sqrt(min(AllMin))<realmin('single')) && ...
                ~(sqrt(max(AllMax))>realmax('single')) && ...
                ~(min(AllMin)^2<realmin('single')) && ...
                ~(max(AllMax)^2>realmax('single'))

        Data.EEG = single(Data.EEG);
		ISr = single(ISr);
        Sources = nan(IShdr.numsolutionpoints,size(Data.EEG,2),3,'single');
    else
    	warning('Temporary conversion to single-precision float might lead to loss of precision in the results !')
    	warning('Trying to keep double-precision float...')
        try
            Sources = nan(IShdr.numsolutionpoints,size(Data.EEG,2),3);
        catch ME
            if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
                warning(ME.message)
            end
            FlagHuge = true;
            warning('Out-of-memory error !')
            warning('Switching back to single-precision floating-point number!')
            warning('Possible loss of precision in the results!')
            Sources = nan(IShdr.numsolutionpoints,size(Data.EEG,2),3,'single');
        end
    end

    if strcmpi(Output,'all')
        % Check if result are likely to produce out-of-memory errors:
        try
            % if ~IShdr.isinversescalar
            Sources = nan(IShdr.numsolutionpoints,size(Data.EEG,2),size(ISr,1)+IShdr.numregularizations);
            % else
        catch ME
            FlagHuge = true;
            if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
                warning(ME.message)
            end
        end
    end

    for reg = 1:IShdr.numregularizations
        fprintf('Regularization %d/%d...\n',reg,IShdr.numregularizations)        
        
        Sources(:,:,1) = squeeze(ISr(1,:,:,reg))*Data.EEG;
        Sources(:,:,2) = squeeze(ISr(2,:,:,reg))*Data.EEG;
        Sources(:,:,3) = squeeze(ISr(3,:,:,reg))*Data.EEG;
        
%         TempFile = fullfile(tempdir,'matlab_compute_inverse',['RISvect_reg',num2str(reg),'.mat']);
%         if exist(fileparts(TempFile),'dir')~=7
%             mkdir(fileparts(TempFile))
%         end
%         save(TempFile,'Sources','-v7.3');
        if strcmpi(Output,'all')
            RISvect{reg} = Sources;
        end
        if ~FlagHuge && strcmpi(Output,'all')
            [RISnorm,RISsvd] = deal(cell(IShdr.numregularizations,1));
        end
        
        if nargout > 1 || strcmpi(Output,'optimal')
            if strcmpi(Output,'all')
                %====== Dipoles norm ======
                fprintf('Computing dipoles norm...\n')
                if ~FlagHuge
                    SourcesNorm = nan(size(Sources,1),size(Sources,2));
                else
                    SourcesNorm = zeros(size(Sources,1),size(Sources,2),'single');
                end
                if size(Sources,2)==1
                    if ~FlagHuge
                        for nsp = 1:size(Sources,1)
                            prc_for_loop(nsp,size(Sources,1),10);
                            SourcesNorm(nsp,:) = sqrt(sum(squeeze(Sources(nsp,:,:)).^2,1));
                        end
                    else
                        for nsp = 1:size(Sources,1)
                            prc_for_loop(nsp,size(Sources,1),10);
                            try
                                SourcesNorm(nsp,:) = sqrt(sum(squeeze(double(Sources(nsp,:,:))).^2,1));
                            catch ME %#ok<NASGU>
                                SourcesNorm(nsp,:) = sqrt(sum(squeeze(Sources(nsp,:,:)).^2,1));
                            end
                        end
                    end
                else
                    if ~FlagHuge
                        for nsp = 1:size(Sources,1)
                            prc_for_loop(nsp,size(Sources,1),10);
                            SourcesNorm(nsp,:) = sqrt(sum(squeeze(Sources(nsp,:,:)).^2,2));
                        end
                    else
                        for nsp = 1:size(Sources,1)
                            prc_for_loop(nsp,size(Sources,1),10);
                            try
                                SourcesNorm(nsp,:) = sqrt(sum(squeeze(double(Sources(nsp,:,:))).^2,2));
                            catch ME %#ok<NASGU>
                                SourcesNorm(nsp,:) = sqrt(sum(squeeze(Sources(nsp,:,:)).^2,2));
                            end
                        end
                    end
                end
            else % strcmpi(Output,'optimal')
                fprintf('Computing average of dipoles norm...\n')
                AvgSourcesNorm = nan(size(Sources,1),1);
                if size(Sources,2)==1
                    for nsp = 1:size(Sources,1)
                        prc_for_loop(nsp,size(Sources,1),10);
                        try
                            AvgSourcesNorm(nsp) = mean(vect(sqrt(sum(squeeze(double(Sources(nsp,:,:))).^2,1))));
                        catch ME %#ok<NASGU>
                            AvgSourcesNorm(nsp) = mean(vect(sqrt(sum(squeeze(Sources(nsp,:,:)).^2,1))));
                        end
                    end
                else
                    for nsp = 1:size(Sources,1)
                        prc_for_loop(nsp,size(Sources,1),10);
                        try
                            AvgSourcesNorm(nsp) = mean(vect(sqrt(sum(squeeze(double(Sources(nsp,:,:))).^2,2))));
                        catch ME %#ok<NASGU>
                            AvgSourcesNorm(nsp) = mean(vect(sqrt(sum(squeeze(Sources(nsp,:,:)).^2,2))));
                        end
                    end
                end
            end
            
            if strcmpi(Output,'optimal')
                fprintf('Computing average norm of dipoles for regularization %d/%d...\n',reg,IShdr.numregularizations)
%                 AvgNorm(reg) = mean(SourcesNorm(:));
                AvgNorm(reg) = mean(AvgSourcesNorm(:));
%             end
%             if strcmpi(Output,'all')
            else % strcmpi(Output,'all')
                if FlagHuge
                    fprintf('Saving temporary solution to disk...\n')
                    TempFile = fullfile(tempdir,'matlab_compute_inverse',['RISnorm_reg',num2str(reg),'.mat']);
                    if exist(fileparts(TempFile),'dir')~=7
                        mkdir(fileparts(TempFile))
                    end
                    save(TempFile,'SourcesNorm','-v7.3');
                    clear SourcesNorm
                else
                    RISnorm{reg} = SourcesNorm;
                end
            end
            
        end
    end
    if strcmpi(Output,'optimal')
        
        TikhonovReg = IShdr.RegularizationValues(:);
        TikhonovReg(TikhonovReg==0)=eps; % otherwise this causes issue with l_corner
        OptimalReg = l_corner(AvgNorm,TikhonovReg);
        
        % Clean memory
        clearvars -except Data ISfile OptimalReg
        
%         % Clean up temporary files
%         rmdir(fileparts(TempFile));
        
        [RISvect,RISnorm,RISsvd] = compute_inverse_no_bad_chans(Data,ISfile,OptimalReg);
        
    else % then Output = 'all'
        if nargout>2
            for reg = 1:IShdr.numregularizations
                fprintf('Regularization %d/%d...\n',reg,IShdr.numregularizations)
                Sources(:,:,1) = squeeze(ISr(1,:,:,reg))*Data.EEG;
                Sources(:,:,2) = squeeze(ISr(2,:,:,reg))*Data.EEG;
                Sources(:,:,3) = squeeze(ISr(3,:,:,reg))*Data.EEG;
                
                %====== First eigenvariate (SVD) ======
                if ~FlagHuge
                    SourcesSVD = nan(size(Sources,1),size(Sources,2));
                else
                    SourcesSVD = zeros(size(Sources,1),size(Sources,2),'single');
                end
                fprintf('Computing SVD of dipoles...\n')
                if size(Sources,2)==1
                    if ~FlagHuge
                        for nsp = 1:size(Sources,1)
                            prc_for_loop(nsp,size(Sources,1),10);
                            Y = squeeze(Sources(nsp,:,:))';
                            
                            % The following is the initial method:
                            % [~,y]=pca(Y,'numcomponents',1);
                            % The following is equivalent, about 10x faster, and matches SPM's eigenvariate
                            [m, n]   = size(Y);
                            [~, s, v] = svd(Y'*Y);
                            s       = diag(s);
                            v       = v(:,1);
                            u       = Y*v/sqrt(s(1));
                            d       = sign(sum(v));
                            u       = u*d;
                            v       = v*d;
                            y       = u*sqrt(s(1)/n);
                            SourcesSVD(nsp,:) = y;
                        end
                    else
                        for nsp = 1:size(Sources,1)
                            prc_for_loop(nsp,size(Sources,1),10);
                            try
                                Y = double(squeeze(Sources(nsp,:,:))');
                            catch ME
                                if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
                                    warning(ME.message)
                                end
                                warning('Running SVD on single-precision sources, possible loss of accuracy...')
                                Y = squeeze(Sources(nsp,:,:))';
                            end
                            
                            % The following is the initial method:
                            % [~,y]=pca(Y,'numcomponents',1);
                            % The following is equivalent, about 10x faster, and matches SPM's eigenvariate
                            [m, n]   = size(Y);
                            [~, s, v] = svd(Y'*Y);
                            s       = diag(s);
                            v       = v(:,1);
                            u       = Y*v/sqrt(s(1));
                            d       = sign(sum(v));
                            u       = u*d;
                            v       = v*d;
                            y       = u*sqrt(s(1)/n);
                            SourcesSVD(nsp,:) = y;
                        end
                    end
                else
                    if ~FlagHuge
                        for nsp = 1:size(Sources,1)
                            prc_for_loop(nsp,size(Sources,1),10);
                            Y = squeeze(Sources(nsp,:,:));
                            
                            % The following is the initial method:
                            % [~,y]=pca(Y,'numcomponents',1);
                            % The following is equivalent, about 10x faster, and matches SPM's eigenvariate
                            [m, n]   = size(Y);
                            [~, s, v] = svd(Y'*Y);
                            s       = diag(s);
                            v       = v(:,1);
                            u       = Y*v/sqrt(s(1));
                            d       = sign(sum(v));
                            u       = u*d;
                            v       = v*d;
                            y       = u*sqrt(s(1)/n);
                            SourcesSVD(nsp,:) = y;
                        end
                    else
                        for nsp = 1:size(Sources,1)
                            prc_for_loop(nsp,size(Sources,1),10);
                            try
                                Y = double(squeeze(Sources(nsp,:,:)));
                            catch ME
                                if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
                                    warning(ME.message)
                                end
                                warning('Running SVD on single-precision sources, possible loss of accuracy...')
                                Y = squeeze(Sources(nsp,:,:));
                            end
                            
                            % The following is the initial method:
                            % [~,y]=pca(Y,'numcomponents',1);
                            % The following is equivalent, about 10x faster, and matches SPM's eigenvariate
                            [m, n]   = size(Y);
                            [~, s, v] = svd(Y'*Y);
                            s       = diag(s);
                            v       = v(:,1);
                            u       = Y*v/sqrt(s(1));
                            d       = sign(sum(v));
                            u       = u*d;
                            v       = v*d;
                            y       = u*sqrt(s(1)/n);
                            SourcesSVD(nsp,:) = y;
                        end
                    end
                end
                RISsvd{reg} = SourcesSVD;
            end
        end
        
%         if FlagHuge
%             fprintf('Loading back norm of dipoles...\n')
%             for reg = 1:IShdr.numregularizations
%                 fprintf('Regularization %d/%d...\n',reg,IShdr.numregularizations)
%                 load(fullfile(tempdir,'matlab_compute_inverse',['RISnorm_reg',num2str(reg),'.mat']),'SourcesNorm')
%                 RISnorm{reg} = SourcesNorm;
%             end
%         end
    end
    
else % sources requested for a particular regularization value only (Output = 'one')
    fprintf('Computing sources in vectorial form...\n')
    
    AllMax = [max(abs(ISr(:))),max(abs(Data.EEG(:)))];
    AllMin = [min(abs(ISr(:))),min(abs(Data.EEG(:)))];
    if ~(sqrt(min(AllMin))<realmin('single')) && ...
                ~(sqrt(max(AllMax))>realmax('single')) && ...
                ~(min(AllMin)^2<realmin('single')) && ...
                ~(max(AllMax)^2>realmax('single'))

        Data.EEG = single(Data.EEG);
		ISr = single(ISr);
        RISvect = nan(IShdr.numsolutionpoints,size(Data.EEG,2),3,'single');
    else
    	warning('Temporary conversion to single-precision float might lead to loss of precision in the results !')
    	warning('Keeping double-precision float... This might lead to out-of-memory errors !')
        try
            RISvect = nan(IShdr.numsolutionpoints,size(Data.EEG,2),3);
        catch ME
            if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
                warning(ME.message)
            end
            FlagHuge = true;
            warning('Out-of-memory error !')
            warning('Switching back to single-precision floating-point number!')
            warning('Likely loss of precision in the results!')
            RISvect = nan(IShdr.numsolutionpoints,size(Data.EEG,2),3,'single');
        end
    end
    
    RISvect(:,:,1) = squeeze(ISr(1,:,:,Regularization))*Data.EEG;
    RISvect(:,:,2) = squeeze(ISr(2,:,:,Regularization))*Data.EEG;
    RISvect(:,:,3) = squeeze(ISr(3,:,:,Regularization))*Data.EEG;
    
    if nargout > 1 % strcmpi(Method,'norm') || strcmpi(Output,'optimal')
        %====== Dipoles norm ======
        fprintf('Computing dipoles norm...\n')
        if ~FlagHuge
            RISnorm = nan(size(RISvect,1),size(RISvect,2));
        else
            RISnorm = zeros(size(RISvect,1),size(RISvect,2),'single');
        end
        if size(RISvect,2)==1
            if ~FlagHuge
                for nsp = 1:size(RISvect,1)
                    prc_for_loop(nsp,size(RISvect,1),10);
                    RISnorm(nsp,:) = sqrt(sum(squeeze(RISvect(nsp,:,:)).^2,1));
                end
            else
                for nsp = 1:size(RISvect,1)
                    prc_for_loop(nsp,size(RISvect,1),10);
                    try
                        RISnorm(nsp,:) = sqrt(sum(squeeze(double(RISvect(nsp,:,:))).^2,1));
                    catch ME %#ok<NASGU>
                        RISnorm(nsp,:) = sqrt(sum(squeeze(RISvect(nsp,:,:)).^2,1));
                    end
                end
            end
        else
            if ~FlagHuge
                for nsp = 1:size(RISvect,1)
                    prc_for_loop(nsp,size(RISvect,1),10);
                    RISnorm(nsp,:) = sqrt(sum(squeeze(RISvect(nsp,:,:)).^2,2));
                end
            else
                for nsp = 1:size(RISvect,1)
                    prc_for_loop(nsp,size(RISvect,1),10);
                    try
                        RISnorm(nsp,:) = sqrt(sum(squeeze(double(RISvect(nsp,:,:))).^2,2));
                    catch ME %#ok<NASGU>
                        RISnorm(nsp,:) = sqrt(sum(squeeze(RISvect(nsp,:,:)).^2,2));
                    end
                end
            end
        end
    end
    if nargout > 2 % elseif strcmpi(Method,'svd')
        if ~FlagHuge
            %====== First eigenvariate (SVD) ======
            RISsvd = nan(size(RISvect,1),size(RISvect,2));
        else
            RISsvd = zeros(size(RISvect,1),size(RISvect,2),'single');
        end
        fprintf('Computing SVD of dipoles...\n')
        if size(RISvect,2)==1
            if ~FlagHuge
                for nsp = 1:size(RISvect,1)
                    prc_for_loop(nsp,size(RISvect,1),10);
                    Y = squeeze(RISvect(nsp,:,:))';
                    
                    % The following is the initial method:
                    % [~,y]=pca(Y,'numcomponents',1);
                    % The following is equivalent, about 10x faster, and matches SPM's eigenvariate
                    [m, n]   = size(Y);
                    [~, s, v] = svd(Y'*Y);
                    s       = diag(s);
                    v       = v(:,1);
                    u       = Y*v/sqrt(s(1));
                    d       = sign(sum(v));
                    u       = u*d;
                    v       = v*d;
                    y       = u*sqrt(s(1)/n);
                    RISsvd(nsp,:) = y;
                end
            else
                for nsp = 1:size(RISvect,1)
                    prc_for_loop(nsp,size(RISvect,1),10);
                    try
                        Y = double(squeeze(RISvect(nsp,:,:))');
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
                            warning(ME.message)
                        end
                        warning('Running SVD on single-precision sources, possible loss of accuracy...')
                        Y = squeeze(RISvect(nsp,:,:))';
                    end
                    
                    % The following is the initial method:
                    % [~,y]=pca(Y,'numcomponents',1);
                    % The following is equivalent, about 10x faster, and matches SPM's eigenvariate
                    [m, n]   = size(Y);
                    [~, s, v] = svd(Y'*Y);
                    s       = diag(s);
                    v       = v(:,1);
                    u       = Y*v/sqrt(s(1));
                    d       = sign(sum(v));
                    u       = u*d;
                    v       = v*d;
                    y       = u*sqrt(s(1)/n);
                    RISsvd(nsp,:) = y;
                end
            end
        else
            if ~FlagHuge
                for nsp = 1:size(RISvect,1)
                    prc_for_loop(nsp,size(RISvect,1),10);
                    Y = squeeze(RISvect(nsp,:,:));
                    % The following is the initial method:
                    % [~,y]=pca(Y,'numcomponents',1);
                    % The following is equivalent, about 10x faster, and matches SPM's eigenvariate
                    [m, n]   = size(Y);
                    [~, s, v] = svd(Y'*Y);
                    s       = diag(s);
                    v       = v(:,1);
                    u       = Y*v/sqrt(s(1));
                    d       = sign(sum(v));
                    u       = u*d;
                    v       = v*d;
                    y       = u*sqrt(s(1)/n);
                    RISsvd(nsp,:) = y;
                end
            else
                for nsp = 1:size(RISvect,1)
                    prc_for_loop(nsp,size(RISvect,1),10);
                    try
                        Y = double(squeeze(RISvect(nsp,:,:)));
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
                            warning(ME.message)
                        end
                        warning('Running SVD on single-precision sources, possible loss of accuracy...')
                        Y = squeeze(RISvect(nsp,:,:));
                    end
                    
                    % The following is the initial method:
                    % [~,y]=pca(Y,'numcomponents',1);
                    % The following is equivalent, about 10x faster, and matches SPM's eigenvariate
                    [m, n]   = size(Y);
                    [~, s, v] = svd(Y'*Y);
                    s       = diag(s);
                    v       = v(:,1);
                    u       = Y*v/sqrt(s(1));
                    d       = sign(sum(v));
                    u       = u*d;
                    v       = v*d;
                    y       = u*sqrt(s(1)/n);
                    RISsvd(nsp,:) = y;
                end
            end
        end
    end
end

end

