
D = loadjson('my_json.json');

for fld = 1:length(Fields)
    for edg = 1:length(D.edge)
        eval(['D.edge{edg}.',Fields{fld},'=[''',char(34),''',num2str(str2double(regexprep(D.edge{edg}.',Fields{fld},',''',char(34),''',''''))+rand*0.1),''',char(34),'''];']);
    end
end

savejson('',D,'my_json_with_rand.json');
