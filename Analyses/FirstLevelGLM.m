function FirstLevelGLM(cfg)
% function FirstLevelGLM(cfg)

nsubjects = length(cfg.subjects);

for sub = 1:nsubjects

    FLdir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.outDir);
    if ~exist(FLdir,'dir'); mkdir(FLdir); end
    scanDir = fullfile(cfg.dataDir,cfg.subjects{sub},cfg.scanDir);
    fprintf('Processing subject %s out of %d \n', cfg.subjects{sub},nsubjects);

    %% Model specification

    if ~exist(fullfile(FLdir,'SPM.mat'),'file')

    % basic info for model set-up
    specification{1}.spm.stats.fmri_spec.dir = {FLdir};
    specification{1}.spm.stats.fmri_spec.timing.units = 'secs';
    specification{1}.spm.stats.fmri_spec.timing.RT = cfg.TR;
    specification{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    specification{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    % load the onsets etc.
    regDir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.regs);
    onsets = str2fullfile(regDir,'run*.mat');
    nSess  = length(onsets);

    % define regressors per session
    for sess = 1:nSess

        fprintf('\t Specifying sessions %d out of %d \n', sess, nSess)
        sessionDir = fullfile(scanDir,sprintf('sess%d',sess));
        regressors = load(onsets{sess});

        % specify the scans
        scans = str2fullfile(sessionDir,cfg.prefix);
        specification{1}.spm.stats.fmri_spec.sess(sess).scans = scans';

        % specify the onsets etc for the different regressors
        if any(cellfun(@isempty,regressors.Onsets)); rm = cellfun(@isempty,regressors.Onsets); regressors.Stim(rm) = []; regressors.Onsets(rm) = []; regressors.Durs(rm) = []; end
        nRegs      = length(regressors.Stim);

        for reg = 1:nRegs
            specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg).name = regressors.Stim{reg};
            specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg).onset = regressors.Onsets{reg};
            specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg).duration = regressors.Durs{reg};
            if contains(cfg.outDir,'viv')
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg).tmod = 0;
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg).pmod.name = 'vividness';
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg).pmod.param = regressors.Viv{reg};
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg).pmod.poly = 1;
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg).orth = 1;
            end
        end

        % specify nuisance regressors
        if ~isempty(regressors.Nuisance)
            if length(regressors.Nuisance)==2; regressors.Nuisance{3} = []; end
            if isempty(regressors.Nuisance{3}); regressors.Nuisance(3) = []; end
            nNuisance = size(regressors.Nuisance,1);
            if isfield(regressors,'Nuisance_names'); nuis_names = regressors.Nuisance_names; else
                nuis_names = {'det_resp','viv_resp','incorrect'}; end
            for n = 1:nNuisance
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg+n).name = nuis_names{n};
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg+n).onset = regressors.Nuisance{n}(:,1);
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg+n).duration = regressors.Nuisance{n}(:,2);
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg+n).tmod = 0;
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg+n).pmod = struct('name', {}, 'param', {}, 'poly', {});
                specification{1}.spm.stats.fmri_spec.sess(sess).cond(reg+n).orth = 1;
            end
        end

        % derivative of realignment parameters
        rpFile = dir(fullfile(sessionDir, 'rp*'));
        fid = fopen(fullfile(rpFile.folder, rpFile.name), 'r');
        rp = fscanf(fid, ' %e %e %e %e %e %e', [6, inf])';
        fclose(fid);
        R = rp;
        R(:, 7:12) = [zeros(1, 6) ; diff(R)];

        % save movement and physio regressors
        regrFile = fullfile(sessionDir,'reg_file.mat');
        save(regrFile, 'R');

        % add movement and physiological regressors
        specification{1}.spm.stats.fmri_spec.sess(sess).multi = {''};
        specification{1}.spm.stats.fmri_spec.sess(sess).regress = struct('name', {}, 'val', {});
        specification{1}.spm.stats.fmri_spec.sess(sess).multi_reg = {regrFile};

        if contains(cfg.outDir,'Localizer')
            specification{1}.spm.stats.fmri_spec.sess(sess).hpf = 128;
        elseif contains(cfg.outDir,'MT')
            specification{1}.spm.stats.fmri_spec.sess(sess).hpf = 256; % lower HPF for task
        else
            error('WRONG TYPE \n')
        end

        % add WM and CSF regressors
        if cfg.prefix(1)=='s'
            [~,wm_mask] = read_nii(fullfile(cfg.dataDir,'rwm_mask.nii'));
            [~,csf_mask] = read_nii(fullfile(cfg.dataDir,'rcsf_mask.nii'));
        elseif cfg.prefix(1)=='u'
            wm_file = str2fullfile(fullfile(cfg.dataDir,cfg.subjects{sub},'Structural'),'c2*.nii');
            csf_file = str2fullfile(fullfile(cfg.dataDir,cfg.subjects{sub},'Structural'),'c3*.nii');

            rwm_file = str2fullfile(fullfile(cfg.dataDir,cfg.subjects{sub},'Structural'),'rc2*.nii');
            if isempty(rwm_file) % reslice
                clear matlabbatch
                matlabbatch{1}.spm.spatial.coreg.write.ref = scans(1);
                matlabbatch{1}.spm.spatial.coreg.write.source = {wm_file;csf_file};
                matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

                addpath(genpath(cfg.spm_dir))
                spm_jobman('run',matlabbatch)
                rmpath(genpath(cfg.spm_dir))
                addpath(cfg.spm_dir)
            end

            rwm_file = str2fullfile(fullfile(cfg.dataDir,cfg.subjects{sub},'Structural'),'rc2*.nii');
            rcsf_file =  str2fullfile(fullfile(cfg.dataDir,cfg.subjects{sub},'Structural'),'rc3*.nii');

            [~,wm_mask] = read_nii(rwm_file);
            [~,csf_mask] = read_nii(rcsf_file);
        end

        wm = zeros(length(scans),1); csf = zeros(length(scans),1);
        for n = 1:length(scans)
            if mod(n,10) == 0
                fprintf('Calculating wm and csf for scan %d out of %d \n',n,length(scans))
            end
            [~,scan] = read_nii(scans{n});
            wm(n)    = mean(scan(wm_mask>0.7));
            csf(n)   = mean(scan(csf_mask>0.7));
            clear scan
        end
        specification{1}.spm.stats.fmri_spec.sess(sess).regress(1).name = 'wm';
        specification{1}.spm.stats.fmri_spec.sess(sess).regress(1).val = wm;
        specification{1}.spm.stats.fmri_spec.sess(sess).regress(2).name = 'csf';
        specification{1}.spm.stats.fmri_spec.sess(sess).regress(2).val = csf;
    end

    specification{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    specification{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0]; % add first derivative
    specification{1}.spm.stats.fmri_spec.volt = 1;
    specification{1}.spm.stats.fmri_spec.global = 'None';
    specification{1}.spm.stats.fmri_spec.mthresh = 0.8;
    if contains(cfg.prefix,'w')  % normalized mask
        specification{1}.spm.stats.fmri_spec.mask = {'D:\NIMADET\Data\gm_mask_logical.nii'};
    else % subject specific mask
        specification{1}.spm.stats.fmri_spec.mask = {''};
    end
    specification{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    % run the model
    addpath(genpath(cfg.spm_dir))
    spm_jobman('run',specification)
    rmpath(genpath(cfg.spm_dir))
    addpath(cfg.spm_dir)
    end

    if ~exist(fullfile(FLdir,'beta_0001.nii'),'file')

    % estimate the model
    spm_file = str2fullfile(FLdir,'SPM.mat');
    estimation{1}.spm.stats.fmri_est.spmmat = {spm_file};
    estimation{1}.spm.stats.fmri_est.method.Classical = 1;

    addpath(genpath(cfg.spm_dir))
    spm_jobman('run',estimation)
    rmpath(genpath(cfg.spm_dir))
    addpath(cfg.spm_dir)
    end
end