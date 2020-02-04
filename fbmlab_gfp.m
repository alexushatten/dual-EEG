%--------------------------------------------------------------------------
% This function calculate the Global Field Power (GFP)
%
% Input : EEG
% 
% Output : GFP
%
%--------------------------------------------------------------------------
% refacto by Renaud Marquis @ FBM lab, October 2018

function GFP=fbmlab_gfp(EEG)

% Calculate GFP
if size(EEG,2)<500000
    GFP=std(EEG);
else
    % Do sections to save memory
%     GFP=[];  % #RM@FBMlab: see below, pre-allocating for speed
    tt=0:1/10:1;
    tt=fix(tt.*size(EEG,2));
    GFP=zeros(sum(diff(tt)),1); % #RM@FBMlab: pre-allocating for speed
    for i=1:numel(tt)-1
%         GFP=[GFP std(EEG(:,tt(i)+1:tt(i+1)))];
        GFP(tt(i)+1:tt(i+1)) = std(EEG(:,tt(i)+1:tt(i+1)));  % #RM@FBMlab: pre-allocating for speed
    end
    GFP = GFP';  % #RM@FBMlab: pre-allocating for speed
    clearvars tt;
end

return;