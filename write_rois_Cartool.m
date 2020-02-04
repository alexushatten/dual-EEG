function ROIs = write_rois_Cartool( Filename, ROIs )
% write_rois_Cartool: write Cartool rois file
%
% write_rois_Cartool( Filename, ROIs )
%
%  Inputs
% --------
% Filename: file path to .rois file to be created
%
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

fileID = fopen(Filename,'w');
fprintf(fileID,'%s\r\n',ROIs.Header);
fprintf(fileID,'%s\r\n',num2str(ROIs.DimensionOrigData));
fprintf(fileID,'%s\r\n',num2str(ROIs.NumberOfROIs));
for n = 1:ROIs.NumberOfROIs
    fprintf(fileID,'%s\r\n',ROIs.ROIname{n});
    fprintf(fileID,'%s\r\n',num2str(ROIs.SolPoints{n}'));
end
fclose(fileID);

end
