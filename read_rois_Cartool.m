function ROIs = read_rois_Cartool( Filename )
% read_roi_Cartool: read Cartool's rois file
%
% ROIs = read_roi_Cartool( Filename )
%
%  Inputs
% --------
% Filename: file path to .rois file
%
%  Outputs
% ---------
% ROIs: structure with :
%                Header: string, magic number, should in principle be 'RO01'
%     DimensionOrigData: double indicating the number of electrodes or
%                        solution points of the original data the ROIs
%                        were made from
%          NumberOfROIs: number of ROIs
%               ROIname: cell array of strings, names of ROIs
%             SolPoints: cell array with vector indicating the solution
%                        points in each ROI
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2019
%-------------------------------------------------------------------------

fid = fopen(Filename);
Count = 1;
while (feof(fid)==0)    % while  not at the end of the file
    line = fgetl(fid); %take line by line
    if Count>3
        if mod(Count,2)==0
            ROIname{Count-3} = line; %#ok<*AGROW,*NASGU>
        else
            TempSolPoints = regexp(line,' ','split');
            SolPoints{Count-4} = str2double(TempSolPoints(~cellfun(@isempty,TempSolPoints))');
        end
    else
        if Count == 1
            if ~strcmp(line,'RO01')
                warning('Magic number ''RO01'' does not match in header')
            end
            Header = line;
        elseif Count == 2
            DimensionOrigData = str2double(line);
        elseif Count == 3
            NumberOfROIs = str2double(line);
        end
    end
    Count = Count+1;
end
fclose(fid);

ROIs.Header = Header;
ROIs.DimensionOrigData = DimensionOrigData;
ROIs.NumberOfROIs = NumberOfROIs;
ROIs.ROIname = ROIname(1:2:end)';
ROIs.SolPoints = SolPoints(1:2:end)';

end
