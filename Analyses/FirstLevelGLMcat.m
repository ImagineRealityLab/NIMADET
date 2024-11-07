function FirstLevelGLMcat(cfg)
% function FirstLevelGLMcat(cfg)
% this one concatenates all runs so that you have one regressor for the
% whole experiment. This is especially useful for things like vividness
% effects which might not have enough variation per run

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
        load(onsets{1},'Stim','Nuisance_names');
        nRegs  = length(Stim);

        % concatenate the onsets
        nScans = nan(nSess,1); Onsets = cell(nRegs,1); Durs = cell(nRegs,1); RJ = cell(nRegs,1);
        Nuisance = cell(length(Nuisance_names),1); Viv = cell(nRegs,1); scans = [];
        for sess = 1:nSess

            sessionDir = fullfile(scanDir,sprintf('sess%d',sess));
            regressors = load(onsets{sess});
            
            tmp = str2fullfile(sessionDir,cfg.prefix);
            nScans(sess) = length(tmp); % count number of scans per session
            scans = cat(2,scans,tmp); clear tmp

            for r = 1:nRegs
                Onsets{r} = cat(1,Onsets{r},regressors.Onsets{r}+(cfg.TR*sum(nScans(1:sess-1))));
                Durs{r} = cat(1,Durs{r},regressors.Durs{r});
                Viv{r}  = cat(1,Viv{r},regressors.Viv{r});
                RJ{r}   = cat(1,RJ{r},regressors.RJ{r});
            end

            for n = 1:length(Nuisance_names)
                if sess == 1; Nuisance{n} = regressors.Nuisance{n};
                else
                    tmp2 = cat(1,Nuisance{n}(:,2),regressors.Nuisance{n}(:,2));
                    tmp1 = cat(1,Nuisance{n}(:,1),regressors.Nuisance{n}(:,1)+(cfg.TR*sum(nScans(1:sess-1))));
                    Nuisance{n} = [tmp1,tmp2];
                end
            end
        end
        
        % define regressors in SPM
        specification{1}.spm.stats.fmri_spec.sess(1).scans = scans';

        % specify the onsets etc for the different regressors
        if any(cellfun(@isempty,Onsets)); rm = cellfun(@isempty,Onsets); Stim(rm) = []; Onsets(rm) = []; Durs(rm) = []; end
        
        for reg = 1:nRegs
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).name = Stim{reg};
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).onset = Onsets{reg};
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).duration = Durs{reg};
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).tmod = 0;

            %specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).pmod(1) = struct('name', {}, 'param', {}, 'poly', {});
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).pmod(1).name = 'vividness';
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).pmod(1).param = Viv{reg};
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).pmod(1).poly = 1;

            %specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).pmod(2) = struct('name', {}, 'param', {}, 'poly', {});
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).pmod(2).name = 'realityjudgement';
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).pmod(2).param = RJ{reg};
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).pmod(2).poly = 1;

            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg).orth = 0; % don't orthogonalize
        end

        % specify nuisance regressors
        if isempty(Nuisance{3}); Nuisance(3) = []; end
        nNuisance = size(Nuisance,1);
        if isfield(regressors,'Nuisance_names'); nuis_names = Nuisance_names; else
            nuis_names = {'det_resp','viv_resp','incorrect'}; end
        for n = 1:nNuisance
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg+n).name = nuis_names{n};
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg+n).onset = Nuisance{n}(:,1);
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg+n).duration = Nuisance{n}(:,2);
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg+n).tmod = 0;
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg+n).pmod = struct('name', {}, 'param', {}, 'poly', {});
            specification{1}.spm.stats.fmri_spec.sess(1).cond(reg+n).orth = 1;
        end

        % create derivative of realignment parameters for all runs
        % concatenated
        R = [];
        for sess = 1:nSess
            sessionDir = fullfile(scanDir,sprintf('sess%d',sess));           
            rpFile = dir(fullfile(sessionDir, 'rp*'));
            fid = fopen(fullfile(rpFile.folder, rpFile.name), 'r');
            rp = fscanf(fid, ' %e %e %e %e %e %e', [6, inf])';
            fclose(fid);
            rp(:, 7:12) = [zeros(1, 6) ; diff(rp)];
            R = cat(1,R,rp);
        end

        % save movement and physio regressors
        regrFile = fullfile(FLdir,'reg_file.mat');
        save(regrFile, 'R');

        % add movement and physiological regressors
        specification{1}.spm.stats.fmri_spec.sess(1).multi = {''};
        specification{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
        specification{1}.spm.stats.fmri_spec.sess(1).multi_reg = {regrFile};
        specification{1}.spm.stats.fmri_spec.sess(1).hpf = 256; % LOWER HPF!

        % add WM and CSF regressors
        [~,wm_mask] = read_nii(fullfile(cfg.dataDir,'rwm_mask.nii'));
        [~,csf_mask] = read_nii(fullfile(cfg.dataDir,'rcsf_mask.nii'));
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
        specification{1}.spm.stats.fmri_spec.sess(1).regress(1).name = 'wm';
        specification{1}.spm.stats.fmri_spec.sess(1).regress(1).val = wm;
        specification{1}.spm.stats.fmri_spec.sess(1).regress(2).name = 'csf';
        specification{1}.spm.stats.fmri_spec.sess(1).regress(2).val = csf;

        % define mask
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
    
    % concatenate the sessions
    addpath(genpath(cfg.spm_dir))
    spm_fmri_concatenate(fullfile(FLdir,'SPM.mat'), nScans');
    rmpath(genpath(cfg.spm_dir))

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