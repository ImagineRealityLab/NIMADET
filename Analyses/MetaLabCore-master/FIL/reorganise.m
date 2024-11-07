function [  ] = reorganise(old_path, new_path, pattern, remove_old,n_dum)
% reorganise(old_path, new_path, old_files, new_files)
% Written by Matan Mazor 2018 
% Adapted by Oliver Warrington 2020

if exist(old_path,'dir')==7
    fname     = spm_select('List',old_path, pattern);
    old_files = cellfun(@(path) fullfile(old_path,path),cellstr(fname),'UniformOutput',0);
    new_files = cellfun(@(path) fullfile(new_path,path),cellstr(fname),'UniformOutput',0);
    mkdir(new_path);
    for i = 1:size(old_files,1)
        copyfile(old_files{i},new_files{i});
    end
    if remove_old; rmdir(old_path,'s'); end

    if (contains(new_path, 'Functional') || contains(new_path, 'Func_loc')) && ~contains(new_path, 'Fieldmaps')
        cd(new_path)
        for d = 1:n_dum
            delete(fname(d, :));
            fprintf('Deleted dummy scan %s.\n',fname(d,:));
        end        
    end
end

end