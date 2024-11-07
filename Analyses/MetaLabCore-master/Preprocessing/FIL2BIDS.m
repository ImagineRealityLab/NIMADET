function [  ] = FIL2BIDS (source_dir, target_dir, sequence_names, sequence_idxs,subj_id)
%FIL2BIDS Gets as input a list of directories and rearranges them in BIDS
%format.
%input: 
%-source_dir: the main subject directory, as retrieved from Charm
%-target_dir: the subject directory, usually named sub-[subject_identifier]
%-sequence_names: a cell array with the sequence types to convert. for
%example, {'anat','func','fmap'};
%=sequence_idxs: a cell array with the corresponding sequence indices. For
%example, {2, 6:10, 3:4};


%Matan Mazor, 2018

%dependencies:
addpath('D:\dicm2nii-master');
setpref('dicm2nii_gui_para', 'save_json', true);

for i=1:length(sequence_names)
    
    sequence_name = sequence_names{i};
    sequence_idx = sequence_idxs{i};
    
    %% 0. Sanity checks
    if sum(strcmp(sequence_name, {'func','fmap','anat'}))==0
        error('Unknown sequence name')
    end

    %% 1. Create subdirectories
    if ~exist(fullfile(target_dir, sequence_name))
        mkdir(fullfile(target_dir, sequence_name))
    end

    %% 2. Untar directories
    mkdir(fullfile(source_dir, sequence_name));
    for i=1:length(sequence_idx)
        untar(fullfile(source_dir, strcat(subj_id,'_FIL.S',num2str(sequence_idx(i)),'.tar')),...
            fullfile(source_dir, sequence_name));
    end

    %% 3. Save as niftis and jsons
    dicm2nii(fullfile(source_dir, sequence_name),fullfile(target_dir, sequence_name),'nii')
    rmdir(fullfile(source_dir, sequence_name),'s');

    %% 4. Validate that everything makes sense, and then change file names
    switch sequence_name
        case 'anat'
            %here I assume that the anatomical image should be mprage
            if ~exist(fullfile(target_dir, sequence_name, 'MPRAGE_64ch_Head.json'), 'file')
                error('Couldn''t find MPRAGE_64ch_Head.json')
            elseif ~exist(fullfile(target_dir, sequence_name, 'MPRAGE_64ch_Head.nii'), 'file')
                error('Couldn''t find MPRAGE_64ch_Head.nii.gz')
            end

            movefile(fullfile(target_dir, sequence_name, 'MPRAGE_64ch_Head.json'),...
                fullfile(target_dir, sequence_name, strcat('sub-',subj_id,'_T1w.json')));

            movefile(fullfile(target_dir, sequence_name, 'MPRAGE_64ch_Head.nii'),...
                fullfile(target_dir, sequence_name, strcat('sub-',subj_id,'_T1w.nii')));

        case 'fmap'
            %here I assume that the fieldmap image is of a certain type
            if ~exist(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_e1.json'), 'file')
                error('Couldn''t find gre_field_mapping_1acq_rl_e1.json')
            elseif ~exist(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_e1.nii'), 'file')
                error('Couldn''t find gre_field_mapping_1acq_rl_e1.nii')
            elseif ~exist(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_e2.json'), 'file')
                error('Couldn''t find gre_field_mapping_1acq_rl_e2.json')
            elseif ~exist(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_e2.nii'), 'file')
                error('Couldn''t find gre_field_mapping_1acq_rl_e2.nii')
            elseif ~exist(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_phase.json'), 'file')
                error('Couldn''t find gre_field_mapping_1acq_rl_phase.json')
            elseif ~exist(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_phase.nii'), 'file')
                error('Couldn''t find gre_field_mapping_1acq_rl_phase.nii')
            end

            movefile(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_e1.json'),...
                fullfile(target_dir, sequence_name, strcat('sub-',subj_id,'_magnitude1_e1.json')));

            movefile(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_e1.nii'),...
                fullfile(target_dir, sequence_name, strcat('sub-',subj_id,'_magnitude1_e1.nii')));
            
            movefile(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_e2.json'),...
                fullfile(target_dir, sequence_name, strcat('sub-',subj_id,'_magnitude1_e2.json')));

            movefile(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_e2.nii'),...
                fullfile(target_dir, sequence_name, strcat('sub-',subj_id,'_magnitude1_e2.nii')));
            
            movefile(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_phase.json'),...
                fullfile(target_dir, sequence_name, strcat('sub-',subj_id,'_phasediff.json')));

            movefile(fullfile(target_dir, sequence_name, 'gre_field_mapping_1acq_rl_phase.nii'),...
                fullfile(target_dir, sequence_name, strcat('sub-',subj_id,'_phasediff.nii')));
        case 'func'
            %here I assume that the _events.tsv files have already been placed in the
            %folder and that they are ordered like the functional scannings
            tsv_files = dir(fullfile(target_dir,'func','*_events.tsv'));
            json_files = dir(fullfile(target_dir,'func','*.json'));
            image_files = dir(fullfile(target_dir,'func','*.nii'));
            if size(tsv_files) ~= size(json_files) 
                error('different number of event files and functional scans')
            end

            for i=1:length(tsv_files)
                run_name = tsv_files(i).name(1:regexp(tsv_files(i).name,'_events.tsv')-1);
                % change json file name
                movefile(fullfile(target_dir, sequence_name, json_files(i).name),...
                           fullfile(target_dir, sequence_name, strcat(run_name,'_bold.json'))); 
                movefile(fullfile(target_dir, sequence_name, image_files(i).name),...
                           fullfile(target_dir, sequence_name, strcat(run_name,'_bold.nii'))); 
            end

    end

end