function open_in_notepadpp(Filename)

if ~ischar(Filename)
    warning('Please input a string')
else
    Path2notepadpp = 'path_to_Notepad++_exe';
    if exist(Path2notepadpp,'file')~=2
        [NotepadppFile,Path2notepadpp] = uigetfile('*.exe','Please select Notepad++ executable');
    else
        [Path2notepadpp,NotepadppFile] = fileparts(Path2notepadpp);
    end
    if NotepadppFile~=0
        system([fullfile(Path2notepadpp,NotepadppFile),' "',Filename,'"']);
    end
end

end