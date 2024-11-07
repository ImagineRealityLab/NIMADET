function fil_mri_organise_from_BIDS(cfg,subj)
% function fil_mri_organise(cfg)
% Sorts folders and deletes dummy volumes
% Details for each subject must be entered in 'fil_subject_details'
% If the data for a subject has already been sorted, then that subject
% is skipped
%
% Receives a cfg struct with required fields:
% cfg.which_subjects = [1:10];
% cfg.which_scans = {'struc', 'func', 'fm',};
% cfg.dir_root = fullfile(path, to, project);
% cfg.n_fun = 6;
% cfg.n_dum = 5;
% cfg.spm_dir = fullfile(path, to, spm12);
%
% Steve Fleming & Dan Bang, FIL, 07/06/2016
% Adapted by Oliver Warrington, FIL, 2020
% Adapted by Nadine Dijkstra to read BIDS format data, FIL, Sept. 2021

%% Directory paths and targets
cwd = pwd;
dir_func    = 'Functional';
dir_block   = 'Run';
dir_struct  = 'Structural';
dir_fm      = 'Fieldmaps';
remove_old  = 0;

%% Reorganise files
for i_s = cfg.which_subjects
    BIDS_root = fullfile(cfg.dir_root, subj{i_s}.name,'MRI','BIDS');
    new_root = fullfile(cfg.dir_root, subj{i_s}.name,'MRI','FIL');
    if ~exist(new_root,'dir'); mkdir(new_root); end

    for scan = cfg.which_scans
        switch scan{:}
            case 'struc'
                old_path = str2fullfile(fullfile(BIDS_root, 'anat'),'*.nii');
                new_path = fullfile(new_root, dir_struct);
                pattern = '^s.*\.nii$';
                reorganise(old_path, new_path, pattern, remove_old);
            case 'func'
                pattern = '^f.*\.nii$';
                for j = 1:length(subj{i_s}.functional)
                    old_path = [old_root, num2str(subj{i_s}.functional(j))];
                    new_path = fullfile(new_root, dir_func, [dir_block, num2str(j)]);
                    reorganise(old_path, new_path, pattern, remove_old,cfg.n_dum);                    
                end                
            case 'fl' % do same for functional localizer as for functional
                pattern = '^f.*\.nii$';
                for j = 1:length(subj{i_s}.func_loc)
                    old_path = [old_root, num2str(subj{i_s}.func_loc(j))];
                    new_path = fullfile(new_root, dir_fl, [dir_block, num2str(j)]);
                    reorganise(old_path, new_path, pattern, remove_old,cfg.n_dum);                    
                end                
            case 'fm'
                pattern = '^s.*\.nii$';
                last = length(subj{i_s}.functional);
                
                for j = length(subj{i_s}.fieldmaps):-2:2
                    for k = length(subj{i_s}.functional(1:last)):-1:1
                        
                        if subj{i_s}.functional(k) > subj{i_s}.fieldmaps(j)
                            last = k-1;
                            
                            % Maps 1 & 2 (phase)
                            old_path = [old_root, num2str(subj{i_s}.fieldmaps(j-1))];
                            new_path = fullfile(new_root, dir_fm, dir_func, [dir_block, num2str(k)]);
                            reorganise(old_path, new_path, pattern, remove_old);
                            
                            % Map 3 (magnitude)
                            old_path = [old_root, num2str(subj{i_s}.fieldmaps(j))];
                            new_path = fullfile(new_root, dir_fm, dir_func, [dir_block, num2str(k)]);
                            reorganise(old_path, new_path, pattern, remove_old);
                        end
                    end
                end
            otherwise
                fprintf('Unrecognised scan. Please check expected "cfg.which_scans" in this script');
        end
        cd(cwd);
    end
end

end