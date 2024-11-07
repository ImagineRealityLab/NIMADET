function convert_4D(which_subjects, root, subj)
% Make 4D niftis to view a movie in FSLeyes

addpath('D:\spm12\matlabbatch');

for SID = which_subjects
    
    fprintf('%s - CONVERTING TO 4D\n', subj{SID}.name);
    nRuns = length(subj{SID}.functional);
    for run = 1:nRuns
        
        fprintf('RUN %d:\n', run);
        
        % Find functional scans
        funcFolder = fullfile(root, subj{SID}.name,'MRI','Functional', sprintf('Run%d', run));
        funcFiles = dir(funcFolder);
        cd(funcFolder);
        
        % Select the 3D niftis
        niiFiles = spm_select('List', fullfile(funcFiles(1).folder, funcFiles(1).name, 'f*'));
        
        % Create batch job
        matlabbatch{1}.spm.util.cat.vols = cellstr(niiFiles);
        matlabbatch{1}.spm.util.cat.name = sprintf('4D_sess%d.nii', run);
        matlabbatch{1}.spm.util.cat.dtype = 0;
        matlabbatch{1}.spm.util.cat.RT = 3.36;
        
        % Run the job
        spm_jobman('initcfg');
        spm_jobman('run', matlabbatch);
    end
end

cd(root)
    