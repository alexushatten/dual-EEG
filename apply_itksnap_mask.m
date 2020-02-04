function apply_itksnap_mask( MaskImagePath, Image2maskPath )
% Remove VOI from image given manually delineated VOI
% to mask made using ITK-SNAP
%
% Usage:
%-------------------------------------------------------------------------
%
% apply_itksnap_mask( MaskImagePath, Image2maskPath )
%
%
% Inputs:
%-------------------------------------------------------------------------
%
% MaskImagePath: path to mask image
%
% Image2maskPath: path to image to msak
%
% Outputs:
%-------------------------------------------------------------------------
%
% Image2maskPath sufixed with "_masked" and where VOI has been masked
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2018
%-------------------------------------------------------------------------

Vtemp = spm_vol(MaskImagePath);

if ~isempty(regexp(Vtemp.descrip,'Cartool','once'))
    itksnap2cartool( MaskImagePath );
    Vmask = spm_vol(spm_file(MaskImagePath,'prefix','r'));
    FlagDelete = true;
else
    Vmask = spm_vol(MaskImagePath);
    FlagDelete = false;
end

Dmask = spm_read_vols(Vmask);

V = spm_vol(Image2maskPath);
D = spm_read_vols(V);

Vnew = V;
Vnew.fname = spm_file(Vnew.fname, 'suffix', '_masked');

Dnew = D.*~Dmask;

spm_write_vol(Vnew, Dnew);

if FlagDelete
    delete(Vmask.fname);
    [p,n,e] = fileparts(Vmask.fname);
    if strcmp(e,'.img') && exist([p,filesep,n,'.hdr'],'file')==2
        if isempty(p)
            delete([n,'.hdr'])
        else
            delete([p,filesep,n,'.hdr']);
        end
    end
end
% flip_mri_left_right( Vnew.fname );

end