function f = normalize_data(seeg)
    f = seeg;
    f.data = [];
    for i=1:length(seeg.header.electrodes_lab)
        signal = seeg.data(i,:);
        normalized = normalize(signal);
    end
    f = normalized
end