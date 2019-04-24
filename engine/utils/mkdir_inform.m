%% creates folder if it does not exist and informs user what happened
function mkdir_inform(nFolder)

if exist(nFolder, 'dir')
    fprintf(['Will write to existing folder: ' nFolder '\n']);
else
    mkdir(nFolder);
    fprintf(['Created folder: ' nFolder '\n']);
end

end