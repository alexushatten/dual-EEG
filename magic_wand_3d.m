function BW3d = magic_wand_3d( ImagePath, OutputName, Overwrite, Tol, MaxDist, FillHoles, Complexity, BiasCorrect )
%magic_wand_3d: perform region growing on 3D image
%
% Usage:
%-------------------------------------------------------------------------
% 
% BW3d = magic_wand_3d( ImagePath, OutputName, Overwrite, Tol, MaxDist,...
% 						FillHoles, Complexity, BiasCorrect )
%
% Inputs:
%-------------------------------------------------------------------------
%
% ImagePath:  full file path to MRI
%
% OutputName: filename for output image, default = 'regiongrowing3d_output'
%
% Overwrite:  if an output of this function is found, whether to replace
%             ('true') or combine ('false') the two runs, default = false.
%
% Tol:        initial tolerance value, e.g 0.1 for 10% of the range across all values in the
%             plane, default = 0.05. After the first detection, though, the
%             tolerance will be re-estimated based on the first shape
%             detected.
%
% MaxDist:    maximal distance from seed point (in pixels), default = Inf
%
% FillHoles: whether to fill holes in the detected clusters, default =
%            false
%
% Complexity: similarity threshold at which the current mask will be
%             extensively examined if it differs from the previous mask,
%             default = 0.5
%
% BiasCorrect: correct for bias field using SPM12 unified segmentation
%              (default = false)
%
% Outputs:
%-------------------------------------------------------------------------
%
% BW3d: Binary mask (with the same size as the input image) indicating
%      1 (true) for associated pixel/voxel and 0 (false) for outside
%
%   In addition, the function writes a binary mask image in the same
%   directory with basename "regiongrowing3d_output".
%
%------------------------------------------------------------------------- 
% NB: SPM required, SPM12 if bias correction requested
% NB: Image processing toolbox required
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2018
%-------------------------------------------------------------------------

%% Init
if nargin < 8
    BiasCorrect = false;
end

if nargin < 7
    Complexity = 0.5;
end

if nargin < 6
    FillHoles = false;
end

if nargin < 5
    MaxDist = Inf;
end

if nargin < 3
    Overwrite = false;
end

if nargin < 2
    OutputName = 'regiongrowing3d_output';
end

%% Check if bias correction requested

if BiasCorrect
    bias_correct(ImagePath);
    ImagePath = spm_file(ImagePath,'prefix','m');
end

%% display image and get seed point

spm_check_registration(ImagePath);
uiwait(msgbox('Click on seed point and when ready click OK'));
xyz = spm_orthviews('Pos',1);
d = spm_get_data(ImagePath,xyz);
V = spm_vol(ImagePath);
% [D,XYZ] = spm_read_vols(V);
D = spm_read_vols(V);

if nargin < 4
    Tol = 0.05*(max(D(:))-min(D(:)));
else
    Tol = Tol*(max(D(:))-min(D(:)));
end

xyz = round(xyz); % for later indexing inside matrix

%% search across planes along 3rd dimension

BW3d = zeros(size(D));
curBW = grayconnected(squeeze(D(:,:,xyz(3))),xyz(1),xyz(2),Tol);
CurSlice = squeeze(D(:,:,xyz(3))); CurSlice(~curBW)=nan;

% Do not trust the user when clicking a representative point because the
% shape to detect is more important than the precise coordinate... So
% recalculate "d" based on what has been detected!

d = prctile(CurSlice(:),50);

% However, do not recalculate Tolerance otherwise there is no way to
% control the sensitivity of the detection...!

% dlim = prctile(CurSlice(:),[20 80]);
% d = prctile(CurSlice(:),50);
% Tol = mean([d-prctile(CurSlice(:),30),prctile(CurSlice(:),70)-d]);

% Just re-calculate the range!
dlim = [d-Tol d+Tol];

% if ~(d >= dlim(1) && d <= dlim(2)),warning('Seed point appears to be not representative of the area detected, results might be unstable'),end
BW3d(:,:,xyz(3)) = curBW;

%% forward
curX = xyz(1);
curY = xyz(2);
curXbkp = curX; curYbkp = curY;
s = xyz(3);
s = s+1;
flagStop = false;
while s <= size(D,3) && ~flagStop
    
    fprintf('Doing slice %d ...\n',s)
    CurSlice = D(:,:,s); CurSlice(~curBW)=nan; % mask data which was not included in previous slice with NaNs
    
    % if any coordinate inside mask of previous slide is in intensity range
    if any((CurSlice(:)>=dlim(1))+(CurSlice(:)<=dlim(2))>1) && ...
        sqrt( (curXbkp-xyz(1))^2 +...
                 (curYbkp-xyz(2))^2 +...
                 (s-xyz(3))^2 ) < MaxDist % and whether we did not reach MaxDist
    
% find coordinate which has both smallest intensity difference and
        % smallest distance to seed point
        [curX, curY] = find((CurSlice>=dlim(1))+(CurSlice<=dlim(2))>1);
        curI = find((CurSlice>=dlim(1))+(CurSlice<=dlim(2))>1);
        Dist = sqrt(sum(([curX, curY]-repmat(xyz(1:2),1,size(curI,1))').^2,2));
        normDist = Dist-min(Dist); normDist = normDist/max(normDist);
        Idiff = abs(CurSlice(curI)-d);
        normIdiff = Idiff-min(Idiff); normIdiff = normIdiff/max(normIdiff);
        [~,theI] = min(sum([normIdiff,normDist],2));
        curX = curX(theI); curY = curY(theI);
        curXbkp = curX; curYbkp = curY;
%    [curX,curY] = find(abs(CurSlice-d)==min(vect(abs(CurSlice-d))),1); % not good, sometimes it goes too far away...

        % backup previous mask
        prevBW = curBW;
        % get mask for new slice by performing 2D region growing based on coordinate
        curBW = grayconnected(D(:,:,s),curX,curY,Tol);
        
        % evalute overlap between current mask and mask of slice s - 1
        Similarity = mean(vect((prevBW - curBW)./(prevBW + curBW)),'omitnan');
        
        if Similarity > Complexity %|| Similarity < -Complexity
            
            % warn user that it's gonna be slow...
            warning('Binary image changed substantially compared to previous iteration... Performing now deep search... This might be relatively slow...')
            
            % if below a certain threshold (related to complexity of structure),
            % , get all coordinates of points that are not in current slice anymore
            % (but were there before in slice s - 1)
            [curX, curY] = find((prevBW - curBW)./(prevBW + curBW)==1);
            
%             [curX, curY] =
%             find((CurSlice>=dlim(1))+(CurSlice<=dlim(2))>1); % too
%             slow...

            curBW = zeros(size(squeeze(D(:,:,s))));
            for x = 1:size(curX,1)
                Diff2d = D(curX(x),curY(x),s)-d;
                if D(curX(x),curY(x),s) >= dlim(1) && D(curX(x),curY(x),s) <= dlim(2)
                    if (Tol-Diff2d)<0
                        curBWtemp = grayconnected(squeeze(D(:,:,s)),curX(x),curY(x),0);
                    else
                        curBWtemp = grayconnected(squeeze(D(:,:,s)),curX(x),curY(x),Tol);
                    end
                    curBW = curBW+curBWtemp;
                end
            end
            curBW = curBW>0;
        end
        BW3d(:,:,s) = curBW;
    else
        curBW = zeros(size(squeeze(D(:,:,s))));
        BW3d(:,:,s) = curBW;
        flagStop = true;
    end
    s = s+1;
end

%% backward

curBW = grayconnected(squeeze(D(:,:,xyz(3))),xyz(1),xyz(2),Tol);
CurSlice = squeeze(D(:,:,xyz(3))); CurSlice(~curBW)=nan;

curX = xyz(1);
curY = xyz(2);
curXbkp = curX; curYbkp = curY;
s = xyz(3);
s = s-1;
flagStop = false;
while s > 0 && ~flagStop
    
    fprintf('Doing slice %d ...\n',s)
    CurSlice = D(:,:,s); CurSlice(~curBW)=nan; % mask data which was not included in previous slice with NaNs
    
    % if any coordinate inside mask of previous slide is in intensity range
    if any((CurSlice(:)>=dlim(1))+(CurSlice(:)<=dlim(2))>1) && ...
        sqrt( (curXbkp-xyz(1))^2 +...
                 (curYbkp-xyz(2))^2 +...
                 (s-xyz(3))^2 ) < MaxDist % and whether we did not reach MaxDist
             
        % find coordinate which has both smallest intensity difference and
        % smallest distance to seed point
        [curX, curY] = find((CurSlice>=dlim(1))+(CurSlice<=dlim(2))>1);
        curI = find((CurSlice>=dlim(1))+(CurSlice<=dlim(2))>1);
        Dist = sqrt(sum(([curX, curY]-repmat(xyz(1:2),1,size(curI,1))').^2,2));
        normDist = Dist-min(Dist); normDist = normDist/max(normDist);
        Idiff = abs(CurSlice(curI)-d);
        normIdiff = Idiff-min(Idiff); normIdiff = normIdiff/max(normIdiff);
        [~,theI] = min(sum([normIdiff,normDist],2));
        curX = curX(theI); curY = curY(theI);
        curXbkp = curX; curYbkp = curY;
%    [curX,curY] = find(abs(CurSlice-d)==min(vect(abs(CurSlice-d))),1); % not good, sometimes it goes too far away...
        
        % backup previous mask
        prevBW = curBW;
        % get mask for new slice by performing 2D region growing based on coordinate
        curBW = grayconnected(D(:,:,s),curX,curY,Tol);
        
        % evalute overlap between current mask and mask of slice s - 1
        Similarity = mean(vect((prevBW - curBW)./(prevBW + curBW)),'omitnan');
        
        if Similarity > Complexity %|| Similarity < -Complexity
        
            % warn user that it's gonna be slow...
            warning('Binary image changed substantially compared to previous iteration... Performing now deep search... This might be relatively slow...')
            
            % if below a certain threshold (related to complexity of structure),
            % , get all coordinates of points that are not in current slice anymore
            % (but were there before in slice s - 1)
            [curX, curY] = find((prevBW - curBW)./(prevBW + curBW)==1);
            
%             [curX, curY] =
%             find((CurSlice>=dlim(1))+(CurSlice<=dlim(2))>1); % too
%             slow...

            curBW = zeros(size(squeeze(D(:,:,s))));
            if size(curX,1)>5000
                uiwait(msgbox('The current slice appears to need more than 5000 iterations of region growing, this might really take a while... Are you sure you want to proceed further?'));
            end
            
            for x = 1:size(curX,1)
                Diff2d = D(curX(x),curY(x),s)-d;
                if D(curX(x),curY(x),s) >= dlim(1) && D(curX(x),curY(x),s) <= dlim(2)
                    if (Tol-Diff2d)<0
                        curBWtemp = grayconnected(squeeze(D(:,:,s)),curX(x),curY(x),0);
                    else
                        curBWtemp = grayconnected(squeeze(D(:,:,s)),curX(x),curY(x),Tol);
                    end
                    curBW = curBW+curBWtemp;
                end
            end
            curBW = curBW>0;
        end
        BW3d(:,:,s) = curBW;
    else
        curBW = zeros(size(squeeze(D(:,:,s))));
        BW3d(:,:,s) = curBW;
        flagStop = true;
    end
    s = s-1;
end

%% fill holes in the masks

if FillHoles
    for s = 1:size(BW3d,3)
        BW3d(:,:,s) = imfill(BW3d(:,:,s),'holes');
    end
end

%% write output image

Vnew = V;
Vnew.fname = spm_file(Vnew.fname,'basename',OutputName);

if exist(Vnew.fname)==2
    if ~Overwrite
        oldBW3d = spm_read_vols(spm_vol(Vnew.fname));
        BW3d = (BW3d + oldBW3d) > 1;
    end
end
spm_write_vol(Vnew,BW3d);

spm_orthviews('AddColouredImage',1, Vnew.fname, [1 0 0]);
% spm_orthviews('AddColouredBlob',1, XYZ, t, V.mat, [1 0 0], 'Region growing output')

end

