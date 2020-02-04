function fix_NifTi4Cartool( ImFilePath, FixImSuffix )
% Fix NifTi image file by flipping / permuting dimensions of the image data
% and corresponding image header so that Cartool can read it and not mess
% it later on...
%
% fix_NifTi4Cartool( ImFilePath, FixImSuffix )
%
%  Inputs
% --------
% ImFilePath  : string, file path to (single) image file to fix (only 1
%               image at a time!)
%
% FixImSuffix : string, suffix to append to filename of input image to make
% the output image
%
%  Outputs
% ---------
% Fixed NifTi / Analyze file converted to Analyze format with suffix
% 'FixImSuffix'.
%
%-------------------------------------------------------------------------
% NB: - Cartool uses the header information to infer orientation flags (e.g.
%     "RAS"), but then "simply to ignore all this mess" (SIC).
%     - We want the exact same image as the input but with an orientation
%     matrix that Cartool is able to "understand" (SIC), so we will need to
%     change both the header and data of the image file, as well as convert
%     it to Analyze format (converting the image file into Analyze via e.g.
%     mricron will "mess" again the image and make it too hard for Cartool
%     to "understand" it).
%     - THIS FUNCTION IS EXPERIMENTAL AND HASN'T GONE THROUGH THOROUGH
%     TESTING YET, USE AT YOUR OWN RISKS!
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, August 2019
%-------------------------------------------------------------------------
if nargin < 2
    FixImSuffix = '';
end

V = spm_vol(ImFilePath);
D = spm_read_vols(V);

% % Get voxel size:
% Scaling = spm_imatrix(V.mat);
% Scaling = Scaling(7:9);

% Find if some axes were "flipped" (in Cartool jargon),
% i.e. "permuted" in Matlab jargon:
[I,Val] = deal(nan(3,1));
for i = 1:3
    [~,I(i)] = max(abs(V.mat(1:3,i)));
    Val(i) = V.mat(I(i),i);
end
% Find if some axes were flipped (aka reversed, see NB in description above):
AxesFlipped = sign(Val)==-1;
% AxesPermuted = I~=J;

%% Fix image header:
% Fix axes:
if isempty(FixImSuffix)
    warning('If ''FixImSuffix'' is empty, this function will overwrite the original image, so I am gonna append ''_fixed'' anyway to save your ass...')
    FixImSuffix = '_fixed';
end
Vnew = V;
Vnew.fname = spm_file(Vnew.fname,'suffix',FixImSuffix,'ext','img');
Vnew.descrip = [Vnew.descrip,' (fixed for Cartool)'];

% Flip (reverse) axes:
Is = find(AxesFlipped);
for i = 1:length(Is)
    P = zeros(1,12);
    P(7:9) = 1;
    P(6+Is(i))=-1;
    S = spm_matrix(P);
    Vnew.mat=Vnew.mat*S;
end

% Permute (flip) axes:
Vnew.mat(1:3,1:3) = Vnew.mat(1:3,I); % should NOT swap rows, just columns, otherwise the transform is not a multiple of eye(3) anymore (because swapping occured twice) !
Vnew.dim = Vnew.dim(I); % because if voxel size is not isotropic spm_write_vol will refuse to write the image

% Roughly (^1) fix origin:
Origin = V.mat(1:3,4);
for i = 1:length(Is)
    Vnew.mat(I(i),4) = sign(Val(i))*V.dim(I(i))+Origin(I(i));
end

% The above works as long as there no rotations, otherwise we need that
% the displacement is not along an axis and cannot be determined solely on
% the base of image dimension and initial origin, but we would need to
% decompose the displacement in terms of sines and cosines for each
% other axes: y and z for x, x and z for y, y and z for x,
% respectively.
% I can get the angles for pitch (x), roll (y) and yaw (z) using:
    % Rot = spm_imatrix(Vnew.mat);
% ... and the rotation parameters would be the 4th to 6th element of the
% Rot variable.
% Though I will eventually implement this in the future but as for now most
% images I need to fix do not have rotations in the transform, I will leave
% it for later.

% (^1) I still could not find the exact transform that fixes the header after image
% data operations (see below), so we will simply align the output image to
% the input image after and eventually find a proper fix.

%% Fix image data:
% Flip (reverse) axes:
for i = 1:length(Is)
    D = flip(D,Is(i));
end
% Permute (flip) axes:
D = permute(D,I);

%% Rewrite image file:
spm_write_vol(Vnew,D);

%% Fix remaining mis-alignement
fprintf('...Refining alignment (fixing approximate origin changes)...\n')
realign_est(Vnew.fname,ImFilePath);

%% Check if rotations are present in the orientation matrix

% In some cases there remain some rotations in the transform matrix,
% but that is usually because there were initially some non-1 or non-0
% values in the initial transform...

end

function realign_est(From,To)

eoptions.quality = 1;
eoptions.sep = 4;
eoptions.fwhm = 1;
eoptions.rtm = 0;
eoptions.interp = 2;
eoptions.wrap = [0 0 0];
eoptions.weight = '';
eoptions.graphics = false;
spm_realign(char({[To,',1'],[From,',1']}),eoptions);

spm_unlink(spm_file(To,'ext','txt','prefix','rp_'));

end

% % From https://en.wikibooks.org/wiki/SPM/How-to#How_to_automatically_reorient_images?
% % Based on https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;d1f675f1.0810
% function auto_reorient(p,To)
% % if ~nargin
% %     [p,sts] = spm_select(Inf,'image');
% %     if ~sts, return; end
% % end
% p = cellstr(p);
% % vg = spm_vol(fullfile(spm('Dir'),'canonical','avg152T1.nii'));
% vg = spm_vol(To);
% tmp = [tempname '.nii'];
% for i=1:numel(p)
%     spm_smooth(p{i},tmp,[12 12 12]);
%     vf = spm_vol(tmp);
%     M  = spm_affreg(vg,vf,struct('regtype','rigid'));
%     [u,s,v] = svd(M(1:3,1:3)); %#ok<ASGLU>
%     M(1:3,1:3) = u*v';
%     N  = nifti(p{i});
%     N.mat = M*N.mat;
%     create(N);
% end
% spm_unlink(tmp);
% end
