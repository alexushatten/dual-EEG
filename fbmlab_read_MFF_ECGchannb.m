function [ECG_channel]=fbmlab_read_MFF_ECGchannb(MFF_filepath)

ECG_channel=0;

Info_XML_file=fullfile(MFF_filepath,'pnsSet.xml');
infoTree=xml_read(Info_XML_file); % from xml_io_tools by Jaroslaw Tuszynski (https://www.mathworks.com/matlabcentral/fileexchange/12907-xml_io_tools)
for i=1:length(infoTree.sensors.sensor)
    if strcmp(infoTree.sensors.sensor(i).name,'ECG')
        ECG_channel=i;
    end   
end

return;