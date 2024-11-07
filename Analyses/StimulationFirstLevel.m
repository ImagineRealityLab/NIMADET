function StimulationFirstLevel(cfg)

nSub = length(cfg.subjects);

addpath(genpath(cfg.spm_dir));

for sub = 1:nSub
    
    fprintf('Processing subject %s \n', cfg.subjects{sub});
    
    % where to save the output
    outDir = fullfile(cfg.results,cfg.subjects{sub},cfg.outDir);
    if ~exist(outDir,'dir'); mkdir(outDir); end
    
    % where the data are
    funcDir = fullfile(cfg.data,cfg.subjects{sub},'Functional');
    regDir  = fullfile(cfg.results,cfg.subjects{sub},'Regressors\Localizer');
    
    % basic info for model set-up
    specification{1}.spm.stats.fmri_spec.dir = {outDir};
    specification{1}.spm.stats.fmri_spec.timing.units = 'secs';
    specification{1}.spm.stats.fmri_spec.timing.RT = cfg.TR;
    specification{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    specification{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
    nSess   = length(cfg.sess);
        
    % define regressors per session
    for sess = 1:nSess
        
        sessionDir = fullfile(funcDir,sprintf('sess%d',cfg.sess(sess)));
        
        % load the onsets
        load(fullfile(regDir,sprintf('run_%d.mat',sess)),...
            'Onsets','Durs');
        
        % specify the scans
        scans = str2fullfile(sessionDir,[cfg.prefix '*.nii']);
        specification{1}.spm.stats.fmri_spec.sess(sess).scans = scans';
        
        % specify the stimulation regressor
        specification{1}.spm.stats.fmri_spec.sess(sess).cond(1).name = 'Stim ON';
        specification{1}.spm.stats.fmri_spec.sess(sess).cond(1).onset = Onsets;
        specification{1}.spm.stats.fmri_spec.sess(sess).cond(1).duration = Durs;
        specification{1}.spm.stats.fmri_spec.sess(sess).cond(1).tmod = 0;
        specification{1}.spm.stats.fmri_spec.sess(sess).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
        specification{1}.spm.stats.fmri_spec.sess(sess).cond(1).orth = 1;
        
        % derivative of realignment parameters
        rpFile = dir(fullfile(sessionDir, 'rp*'));
        fid = fopen(fullfile(rpFile.folder, rpFile.name), 'r');
        rp = fscanf(fid, ' %e %e %e %e %e %e', [6, inf])';
        fclose(fid);
        R = rp;
        R(:, 7:12) = [zeros(1, 6) ; diff(R)];
        
        % save movement regressors
        regrFile = fullfile(sessionDir, sprintf('mov_sess%d.mat', sess));
        save(regrFile, 'R');
        fclose('all');
        
        % add movement regressors
        specification{1}.spm.stats.fmri_spec.sess(sess).multi = {''};
        specification{1}.spm.stats.fmri_spec.sess(sess).regress = struct('name', {}, 'val', {});
        specification{1}.spm.stats.fmri_spec.sess(sess).multi_reg = {regrFile};
        specification{1}.spm.stats.fmri_spec.sess(sess).hpf = 128;        
        
        clear scans Onsets Durs
    end
    
    specification{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    specification{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
    specification{1}.spm.stats.fmri_spec.volt = 1;
    specification{1}.spm.stats.fmri_spec.global = 'None';
    specification{1}.spm.stats.fmri_spec.mthresh = 0.8;
    specification{1}.spm.stats.fmri_spec.mask = {''};
    specification{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    % run the model
    spm_jobman('run',specification)
    
    % estimate the model
    spm_file = str2fullfile(outDir,'SPM.mat');
    estimation{1}.spm.stats.fmri_est.spmmat = {spm_file};
    estimation{1}.spm.stats.fmri_est.method.Classical = 1;
    
    spm_jobman('run',estimation)
    
    % get the contrast for stimulus ON
    matlabbatch{1}.spm.stats.con.spmmat = {spm_file};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Stimulation';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 0];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
    matlabbatch{1}.spm.stats.con.delete = 0;
    
    spm_jobman('run', matlabbatch)
end