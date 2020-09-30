function h = remove_specified_data_from_array(list, letter_to_remove)
    h = list;
    h.header.electrodes_lab = [];
    h.data =[];
    h.header.electrodes_index = [];
    array_to_remove = [];
    j = 1;
    for i=1:length(list.header.electrodes_lab)
        if contains(list.header.electrodes_lab(i),letter_to_remove)
            continue
        else
            h.header.electrodes_lab = [h.header.electrodes_lab list.header.electrodes_lab(i)];
            h.data(j,:) = list.data(i,:);
            h.header.electrodes_index = [h.header.electrodes_index j];
            j = j + 1 ;
            
        end
    end
    h.header.electrodes_lab = h.header.electrodes_lab';
    h.header.electrodes_index = h.header.electrodes_index';
end