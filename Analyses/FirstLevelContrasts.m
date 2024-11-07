function FirstLevelContrasts(cfg)

nSubjects = length(cfg.subjects);
if isfield(cfg,'tConNames'); nTCons = length(cfg.tConNames); else; nTCons = 0; end
if isfield(cfg,'fConNames'); nFCons = length(cfg.fConNames); else; nFCons = 0; end

for sub = 1:nSubjects
    
    FLdir = fullfile(cfg.root,cfg.subjects{sub},cfg.FLdir);
    spm_file = fullfile(FLdir,'SPM.mat');    
    
    matlabbatch{1}.spm.stats.con.spmmat = {spm_file};
    load(spm_file,'SPM')
    nCols = length(SPM.xX.name);

    % define the t-contrasts
    if cfg.plotCon; tconMat = zeros(nTCons,nCols); end
    for tcon = 1:nTCons

        tvector = zeros(nCols,1);

        if ~iscellstr(cfg.tCons{tcon})
            tvector(1:length(cfg.tCons{tcon})) = cfg.tCons{tcon};
        else
            if length(cfg.tCons{tcon})==2 && size(cfg.tCons{tcon},1)==1 % one vs other
                cond1 = contains(SPM.xX.name,cfg.tCons{tcon}{1}) & ...
                    contains(SPM.xX.name,'bf(1)');% & ... % first derivative
                    %~contains(SPM.xX.name,'x');

                cond2 = contains(SPM.xX.name,cfg.tCons{tcon}{2}) & ...
                    contains(SPM.xX.name,'bf(1)');% & ... % first derivative
                    %~contains(SPM.xX.name,'x'); % no pmods

                tvector(cond1) = 1;
                tvector(cond2) = -1;        

            elseif length(cfg.tCons{tcon}) == 1 && contains(cfg.tCons{tcon},'x')
    
                c = 1;
                if contains(cfg.tCons{tcon},'stimulus')
                    conds{c} = {'pres','abs'}; c = c+1;
                end

                if contains(cfg.tCons{tcon},'congruency')
                    conds{c} = {'cong','inco'}; c = c+1;
                end

                if contains(cfg.tCons{tcon},'orientation')
                    conds{c} = {'ori1','ori2'}; 
                end

                cond1 = ((contains(SPM.xX.name,conds{1}{1}) &  contains(SPM.xX.name,conds{2}{1})) | ...
                    (contains(SPM.xX.name,conds{1}{2}) &  contains(SPM.xX.name,conds{2}{2}))) & ...
                    contains(SPM.xX.name,'bf(1)') & ... % first derivative
                    ~contains(SPM.xX.name,'x');

                cond2 = ((contains(SPM.xX.name,conds{1}{1}) & contains(SPM.xX.name,conds{2}{2})) | ...
                    (contains(SPM.xX.name,conds{1}{2}) &  contains(SPM.xX.name,conds{2}{1}))) & ...
                    contains(SPM.xX.name,'bf(1)') & ... % first derivative
                    ~contains(SPM.xX.name,'x');


                tvector(cond1) = 1;
                tvector(cond2) = -1;

            elseif contains(cfg.tCons{tcon},'three way') % 3-way interaction

                conds = {'inco_abs_ori1','inco_abs_ori2','inco_pres_ori1','inco_pres_ori2',...
                    'cong_abs_ori1','cong_abs_ori2','cong_pres_ori1','cong_pres_ori2'};

                cong = [-1 -1 -1 -1 1 1 1 1]; stim = [-1 -1 1 1 -1 -1 1 1];
                ori  = [1 -1 1 -1 1 -1 1 -1];
                interaction = cong.*stim.*ori;

                cond1_names = conds(interaction==1);
                cond2_names = conds(interaction==-1);

                cond1_idx = zeros(length(tvector),1);
                for c1 = 1:length(cond1_names)
                    tmp = contains(SPM.xX.name',cond1_names{c1});
                    cond1_idx = cond1_idx+tmp; clear tmp
                end

                cond2_idx = zeros(length(tvector),1);
                for c1 = 1:length(cond2_names)
                    tmp = contains(SPM.xX.name',cond2_names{c1});
                    cond2_idx = cond2_idx+tmp; clear tmp
                end

                cond1 = cond1_idx'==1 &...
                    contains(SPM.xX.name,'bf(1)') &... % first derivative
                    ~contains(SPM.xX.name,'x');

                cond2 = cond2_idx'==1 &...
                    contains(SPM.xX.name,'bf(1)') &... % first derivative
                    ~contains(SPM.xX.name,'x');

                tvector(cond1) = 1;
                tvector(cond2) = -1;


            elseif length(cfg.tCons{tcon}) == 1 % one vs rest
                reg = contains(SPM.xX.name,cfg.tCons{tcon}) & ...
                    contains(SPM.xX.name,'bf(1)');

                tvector(reg) = 1;
            end
        end

        tvector(sum(SPM.xX.X)==0) = 0;


        if cfg.plotCon; tconMat(tcon,:) = tvector; end
        matlabbatch{1}.spm.stats.con.consess{tcon}.tcon.name = cfg.tConNames{tcon};
        matlabbatch{1}.spm.stats.con.consess{tcon}.tcon.weights = tvector;
        matlabbatch{1}.spm.stats.con.consess{tcon}.tcon.sessrep = 'none';
    end

    
    % define the f contrasts
    for fcon = 1:nFCons
        matlabbatch{1}.spm.stats.con.consess{fcon+tcon}.fcon.name = cfg.fConNames{fcon};
        matlabbatch{1}.spm.stats.con.consess{fcon+tcon}.fcon.weights = cfg.fConVec{fcon};
        matlabbatch{1}.spm.stats.con.consess{fcon+tcon}.fcon.sessrep = 'repl';
    end    
    
    % delete if already existing contrasts
    matlabbatch{1}.spm.stats.con.delete = 0; % Do not delete previous contrasts
    
    % plot the contrasts?
    if cfg.plotCon
        figure; 
        imagesc(tconMat);         
    end
    
    addpath(genpath(cfg.spm_dir))
    spm_jobman('run',matlabbatch)
    rmpath(genpath(cfg.spm_dir))
    addpath(cfg.spm_dir);
    
end

