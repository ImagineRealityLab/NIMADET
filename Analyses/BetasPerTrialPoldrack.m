function BetasPerTrialPoldrack(cfg)
% gets betas per trials
% Also added temporal derivatives to take into account the temporal jitter
% between slices caused by the multi-band acquisition which is too small to
% correct with slice-time correction

nSub = length(cfg.subjects);

addpath(genpath(cfg.spm_dir));

for sub = 1:nSub
    
    fprintf('Processing subject %s \n', cfg.subjects{sub});
    
    % where to save the output
    if cfg.prefix(1) == 'w'
    outDir = fullfile(cfg.results,cfg.subjects{sub},'FirstLevel',...
        [cfg.type '_perTrialPoldrack']); if ~exist(outDir,'dir'); mkdir(outDir); end
    elseif cfg.prefix(1) == 's'
    outDir = fullfile(cfg.results,cfg.subjects{sub},'FirstLevel',...
        [cfg.type '_perTrialPoldrack_smooth']); if ~exist(outDir,'dir'); mkdir(outDir); end    
    elseif cfg.prefix(1) == 'u'
        outDir = fullfile(cfg.results,cfg.subjects{sub},'FirstLevel/SS',...
        [cfg.type '_perTrialPoldrack']); if ~exist(outDir,'dir'); mkdir(outDir); end    
    end

    % check if already done
    if exist(fullfile(outDir,'RUN_01'),'dir')
        fprintf('Subject %s already done - moving on \n',cfg.subjects{sub})
    else
    
    % where the data are
    funcDir = fullfile(cfg.data,cfg.subjects{sub},'Functional');
    regDir  = fullfile(cfg.results,cfg.subjects{sub},'Regressors',cfg.type);
    
    % define regressors per session
    nSess   = length(cfg.sess);
    for sess = 1:nSess
        
        fprintf('\t Specifying session %d out of %d \n', sess, nSess)
        
        % load the onsets
        load(fullfile(regDir,sprintf('run_%d.mat',sess)),...
            'Onsets','Durs','Stim','Nuisance');
        if ~exist('Nuisance','var'); Nuisance = []; end
        if sum(isnan(Onsets)); nan_idx = isnan(Onsets); Onsets(nan_idx) = []; Durs(nan_idx) = []; end
        
        % define output dir per run
        sessionDir = fullfile(funcDir,sprintf('sess%d',cfg.sess(sess)));
        sessionoutDir = fullfile(outDir,sprintf('RUN_%02d',sess));
        if ~exist(sessionoutDir,'dir'); mkdir(sessionoutDir); end
        
        % find the scans
        scans = str2fullfile(sessionDir,[cfg.prefix '*.nii']);         
        
        %% Model set-up        
        specification{1}.spm.stats.fmri_spec.timing.units = 'secs';
        specification{1}.spm.stats.fmri_spec.timing.RT = cfg.TR;
        specification{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        specification{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        specification{1}.spm.stats.fmri_spec.sess.scans = scans';
        
        %% Physio, WM + GM regs and bases functions
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
        specification{1}.spm.stats.fmri_spec.sess.multi = {''};
        specification{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        specification{1}.spm.stats.fmri_spec.sess.multi_reg = {regrFile};

        % HPF
        if strcmp(cfg.type,'Localizer')            
            specification{1}.spm.stats.fmri_spec.sess.hpf = 128;
        elseif strcmp(cfg.type,'MT')
            specification{1}.spm.stats.fmri_spec.sess.hpf = 256; % lower HPF for task
        else
            error('WRONG TYPE \n')
        end

        % add WM and CSF regressors
        if cfg.prefix(1)=='w'
            [~,wm_mask] = read_nii(fullfile(cfg.data,'rwm_mask.nii'));
            [~,csf_mask] = read_nii(fullfile(cfg.data,'rcsf_mask.nii'));
        elseif cfg.prefix(1)=='u'
            wm_file = str2fullfile(fullfile(cfg.data,cfg.subjects{sub},'Structural'),'c2*.nii');
            csf_file = str2fullfile(fullfile(cfg.data,cfg.subjects{sub},'Structural'),'c3*.nii');

            rwm_file = str2fullfile(fullfile(cfg.data,cfg.subjects{sub},'Structural'),'rc2*.nii');
            if isempty(rwm_file) % reslice
                clear matlabbatch
                matlabbatch{1}.spm.spatial.coreg.write.ref = scans(1);
                matlabbatch{1}.spm.spatial.coreg.write.source = {wm_file;csf_file};
                matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

                spm_jobman('run',matlabbatch)                
            end

            rwm_file = str2fullfile(fullfile(cfg.data,cfg.subjects{sub},'Structural'),'rc2*.nii');
            rcsf_file =  str2fullfile(fullfile(cfg.data,cfg.subjects{sub},'Structural'),'rc3*.nii');

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
        specification{1}.spm.stats.fmri_spec.sess.regress(1).name = 'wm';
        specification{1}.spm.stats.fmri_spec.sess.regress(1).val = wm;
        specification{1}.spm.stats.fmri_spec.sess.regress(2).name = 'csf';
        specification{1}.spm.stats.fmri_spec.sess.regress(2).val = csf;
        
        % bases functions
        specification{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        specification{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0]; % add 1st derivative 
        specification{1}.spm.stats.fmri_spec.volt = 1;
        specification{1}.spm.stats.fmri_spec.global = 'None';
        specification{1}.spm.stats.fmri_spec.mthresh = 0.8;
        if contains(cfg.prefix,'w')
            specification{1}.spm.stats.fmri_spec.mask = {'D:\NIMADET\Data\gm_mask_logical.nii'};
        else
            specification{1}.spm.stats.fmri_spec.mask = {''};
        end
        specification{1}.spm.stats.fmri_spec.cvi = 'AR(1)';               
        
        %% Run one model per trial
        nTrls = length(Onsets);
        for t = 1:nTrls
            
            fprintf('Estimating trial %d out of %d \n', t,nTrls)
            trl_dir = fullfile(sessionoutDir,sprintf('TRL_%03d',t));
            if ~exist(trl_dir,'dir'); mkdir(trl_dir); end
            
            specification{1}.spm.stats.fmri_spec.dir = {trl_dir};
            
            % regressor for that one trial
            if strcmp(cfg.type,'Localizer')
                specification{1}.spm.stats.fmri_spec.sess.cond(1).name = sprintf('trl_%d_stm_%d',t,Stim(t));
            else
                specification{1}.spm.stats.fmri_spec.sess.cond(1).name = sprintf('trl_%d_%s_',t,Stim{t});
            end
            specification{1}.spm.stats.fmri_spec.sess.cond(1).onset = Onsets(t);
            specification{1}.spm.stats.fmri_spec.sess.cond(1).duration = Durs(t);
            
            % other trials
            trl_idx = setdiff(1:nTrls,t);
            specification{1}.spm.stats.fmri_spec.sess.cond(2).name = 'other_trials';
            specification{1}.spm.stats.fmri_spec.sess.cond(2).onset = Onsets(trl_idx);
            specification{1}.spm.stats.fmri_spec.sess.cond(2).duration = Durs(trl_idx);
            count = 2;
            
            % nuisance regressors
            if ~isempty(Nuisance)
                specification{1}.spm.stats.fmri_spec.sess.cond(count+1).name = 'Detection text';
                specification{1}.spm.stats.fmri_spec.sess.cond(count+1).onset = Nuisance{1}(:,1);
                specification{1}.spm.stats.fmri_spec.sess.cond(count+1).duration = Nuisance{1}(:,2);
                
                specification{1}.spm.stats.fmri_spec.sess.cond(count+2).name = 'Detection response';
                specification{1}.spm.stats.fmri_spec.sess.cond(count+2).onset = Nuisance{1}(:,1)+Nuisance{1}(:,2);
                specification{1}.spm.stats.fmri_spec.sess.cond(count+2).duration = zeros(length(Nuisance{1}),1);
                
                specification{1}.spm.stats.fmri_spec.sess.cond(count+3).name = 'Vividness text';
                specification{1}.spm.stats.fmri_spec.sess.cond(count+3).onset = Nuisance{2}(:,1);
                specification{1}.spm.stats.fmri_spec.sess.cond(count+3).duration = Nuisance{2}(:,2);
                
                specification{1}.spm.stats.fmri_spec.sess.cond(count+4).name = 'Vividness response';
                specification{1}.spm.stats.fmri_spec.sess.cond(count+4).onset = Nuisance{2}(:,1)+Nuisance{2}(:,2);
                specification{1}.spm.stats.fmri_spec.sess.cond(count+4).duration = zeros(length(Nuisance{2}),1);
            end
            
            % run the model
            spm_jobman('run',specification)
            spm_file = str2fullfile(trl_dir,'SPM.mat');
            estimation{1}.spm.stats.fmri_est.spmmat = {spm_file};
            estimation{1}.spm.stats.fmri_est.method.Classical = 1;            
            spm_jobman('run',estimation)
            
            % only keep the trl beta
            trl_beta = fullfile(trl_dir,'beta_0001.nii');
            if strcmp(cfg.type,'Localizer')
                new_name = fullfile(sessionoutDir,sprintf('trl_%d_stm_%d.nii',t,Stim(t)));
            else
                new_name = fullfile(sessionoutDir,sprintf('trl_%d_%s_.nii',t,Stim{t}));
            end
            copyfile(trl_beta,new_name);
            rmdir(trl_dir,'s');
        end
        
        
        %% Clean up after run
        clear scans Onsets Stim Durs specification
    end
    end
end

rmpath(genpath(cfg.spm_dir));
addpath(cfg.spm_dir);