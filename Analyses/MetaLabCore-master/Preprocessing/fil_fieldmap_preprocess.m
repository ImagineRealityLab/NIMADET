function fil_fieldmap_preprocess(which_subjects,subj)
% function fil_fieldmap_preprocess(which_subjects)
%
% Script calling SPM's fieldmap_preprocess for multi-subject, 
% multi-session datasets (NB: reorient images first in SPM).
%
% Fieldmap_preprocess creates a single fieldmap (fpm5_*.), converts it to
% a voxel displacement file (vdm5_*.), and matches the vdm5 file to the
% first EPI image from each run.
%
% For multiple sessions, copies of original fieldmaps should be present in 
% each of dir_fieldmap/session* directories. This allows each session to have
% individual Match VDM to EPI maps.
% 
% See http://intranet.fil.ion.ucl.ac.uk/wiki/physicswiki/doku.php?id=start:data_processing:using_field_maps
% for further information [MUST BE ON FIL NETWORK TO READ]
%
% which_subjects is a vector
%
% Steve Fleming & Dan Bang, FIL, 07/06/2016
% Adapted by Nadine Dijkstra, FIL, Sept. 2021

%% Directory paths and targets
dir_funct   = 'Functional';
dir_block   = 'sess';
dir_fm      = 'Fieldmaps';
dir_root    = 'D:\NIMADET\Data';
dir_spm     = 'D:\spm12';

%% Add SPM directory
addpath(genpath(dir_spm));
rmpath(genpath([dir_spm, '\external\fieldtrip']));

%% Fieldmap_preprocess parameters (see help Fieldmap_preprocess for details) 
% settings for sequence ‘cmrr_mbep2d_bold_R016_2mm_MB4_FoV212’ and the fieldmap ‘gre_field_mapping_1acq_rl’
te1             = 10.0;     % short echo time 
te2             = 12.46;    % long echo time
epifm           = 0;        % epi-based fieldmap (1/0)
tert            = 58.3;       % total echo (EPI) readout time -- see EPI sequence
kdir            = -1;       % blip direction (+1/-1) -- see EPI sequence
mask            = 0;        % (optional flag, default=1) Do brain masking or not (only if non-epi fieldmap)
match           = 1;        % (optional flag, default=1) Match fieldmap to epi or not

% loop through all subjects and sessions
%===========================================================================
for i_s = which_subjects
   
    % display current subject
    fprintf(['====SUBJECT ',num2str(i_s),': fieldmap pre-processing\n']);
    
    % functional
    n_fun = length(subj{i_s}.functional);
    for j = 1:n_fun
        % find folders
        fm_dir  = fullfile(dir_root,subj{i_s}.name,dir_fm,'Functional',[dir_block,num2str(j)]);
        epi_dir = fullfile(dir_root,subj{i_s}.name,dir_funct,[dir_block,num2str(j)]);
        % run fieldmap_preprocess
        FieldMap_preprocess(fm_dir,epi_dir,[te1, te2, epifm, tert, kdir, mask, match]);
        fprintf(['....pre-processing completed for main task block ',num2str(j),'\n']);
    end
    
end

end