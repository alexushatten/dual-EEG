function [thedata,numclusters,numfiles,numtimeframes,numtracks,tracknames] = read_seg_Cartool(openfilename)
% read_seg_Cartool: opens a Cartool segmentation results file (.seg)
%
% [thedata,numclusters,numfiles,numtimeframes,numtracks,tracknames] = read_seg_Cartool(openfilename)
%
%  Inputs
% --------
% openfilename: string, path to .seg file
%
%  Outputs
% ---------
% thedata: [numtracks x numtimeframes] double data array
% numclusters: integer, number of clusters
% numfiles: integer, number of files
% numtimeframes: integer, number of time frames
% numtracks: integer, number of tracks (e.g. 5)
% tracknames: integer, track names (e.g. GFP, Dis, Seg, GEV, Corr)
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018, based on script by pierre.megevand@medecine.unige.ch
%-------------------------------------------------------------------------

% open filename for reading in text mode
fid=fopen(openfilename,'rt');

% for .eph files
if strcmp(openfilename(end-3:end),'.seg')==1
    
    % read header
    eph=textscan(fid,'%s','delimiter','/n');
    eph=eph{1};
    
    numclusters=sscanf(eph{1},'%f',1);
    numfiles=sscanf(eph{1},'%*f %f',1);
    numtimeframes=sscanf(eph{1},'%*f %*f %f',1);
    
    H2 = regexp(eph{2},'\s','split');
    
    numtracks = str2double(H2{1});
    tracknames = H2(2:end);
    
    % prepare for reading data
    formatstring='%f';
    if (numfiles*numtracks)>1
        for i=1:(numfiles*numtracks)-1
            formatstring=[formatstring ' %f'];
        end
    end

    % read data
    thedata=zeros(numtimeframes,(numfiles*numtracks));
    for i=1:numtimeframes
        thedata(i,:)=sscanf(eph{i+2},formatstring);
    end

else
    error('incorrect file type');
end

% close file
fclose(fid);

end