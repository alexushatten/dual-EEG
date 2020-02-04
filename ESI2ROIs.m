function SourcesInROIs = ESI2ROIs( Sources, ROISfilepath )
% ESI2ROIs: perform SVD on ESI to go from solution points to ROIs defined
% by atlas
%
% SourcesInROIs = ESI2ROIs( Sources, Parcellation )

%  Inputs
% --------
% Sources:      [SPs x time] reconstructed sources
% ROISfilepath:  char, path to .rois file (output of Cartool)
%
%  Outputs
% ---------
% SourcesInROIs: [ROIs x time] reconstructed sources using SVD
%
% See also COMPUTE_INVERSE
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2019
%-------------------------------------------------------------------------

ROIs = read_rois_Cartool(ROISfilepath);

try
    SourcesInROIs = nan(ROIs.NumberOfROIs,size(Sources,2));
catch ME
    if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
        warning(ME.message)
    end
    warning('Possible out-of-memory error !')
    warning('Switching to single-precision floating-point number...')
    SourcesInROIs = nan(ROIs.NumberOfROIs,size(Sources,2),'single');
end
for r = 1:ROIs.NumberOfROIs
    try 
        Y = double(Sources(ROIs.SolPoints{r},:))';
    catch ME %#ok<NASGU>
        Y = Sources(ROIs.SolPoints{r},:)';
    end
    [m, n]   = size(Y);
    [~, s, v] = svd(Y'*Y);
    s       = diag(s);
    v       = v(:,1);
    u       = Y*v/sqrt(s(1));
    d       = sign(sum(v));
    u       = u*d;
    v       = v*d;
    y       = u*sqrt(s(1)/n);
    SourcesInROIs(r,:) = y;
end

end