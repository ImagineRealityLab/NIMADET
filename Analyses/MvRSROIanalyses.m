function MvRSROIanalyses(cfg)
% trains and tests LDA classifiers with n-fold cross-validation on the task
% data

rng(1); % for reproduceability
lambda = 100000; % for ridge regularization 

% get ROIs
group_dir = fullfile(cfg.root,'GroupResults',cfg.dir);
nROIs     = length(cfg.ROIs);
for r = 1:nROIs
    [~,ROI] = read_nii(fullfile(group_dir,[cfg.ROIs{r} '.nii']));
    if r ==1
        atlas = ROI; 
    else
        atlas = atlas+(ROI.*r); 
    end
    clear ROI
end
mask = atlas>0; atlas = atlas(mask);

% loop over subjects
nSub = length(cfg.subjects);
acc     = nan(nSub,nROIs,2,2); % per congruency, vividness and RJ
MSE     = nan(nSub,nROIs,2); % just viv and RJ 
for sub = 1:nSub

    fprintf('Processing subject %s \n',  cfg.subjects{sub})

    % task data
    dataDir = fullfile(cfg.data,cfg.subjects{sub},'FirstLevel',['MT' cfg.dataDir]);
    nRuns   = length(dir(dataDir))-2;
    D = []; I = []; P = []; R = []; V = []; data = []; run_idx = [];
    for run = 1:nRuns

        fprintf('Reading in task betas for run %d... \n',run)

        if contains(cfg.dataDir,'Poldrack')
            beta_niftis = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),'trl_*.nii');
            nifti_names = cell(length(beta_niftis),1); % re-sort
            beta_names  = cell(length(nifti_names),1);
            for n = 1:length(beta_niftis)
                nifti_names{n} = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),sprintf('trl_%d_*.nii',n));
                [~,beta_names{n}] = fileparts(nifti_names{n});
            end
        else
            load(fullfile(dataDir,sprintf('RUN_%02d',run),'SPM.mat'),'SPM');
            bidx = find(contains(SPM.xX.name,'bf(1)') & contains(SPM.xX.name,'Ima'));
            beta_niftis = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),'beta_*.nii');
            beta_names = SPM.xX.name(bidx);
            nifti_names = beta_niftis(bidx);
        end

        nTrls = length(beta_names);

        % remove blocks failed imagery check
        if cfg.ima_check
            load(fullfile(fileparts(cfg.data),'Data',cfg.subjects{sub},...
                'Behaviour',sprintf('MT_%s_run%d.mat',cfg.subjects{sub},run)),'C')
            if any(C==0)
                ridx = [];
                for c = 1:length(C)
                    idx = (c-1)*(nTrls/length(C))+1:c*(nTrls/length(C));
                    if C(c)==0
                        ridx = [ridx,idx];
                    end
                end
                beta_names(ridx) = [];
                nifti_names(ridx) = [];
            end
        end

        nTrls = length(beta_names);
        betas = nan(nTrls,sum(mask(:)));
        for t = 1:nTrls
            [hdr2,Y] = read_nii(nifti_names{t});
            betas(t,:) = Y(mask); clear Y;

            num_idx = [strfind(beta_names{t},'Det_') strfind(beta_names{t},'_Ima')];
            D = [D str2double(beta_names{t}(num_idx(1)+4:num_idx(2)-1))];

            num_idx = [strfind(beta_names{t},'Ima_') strfind(beta_names{t},'_P')];
            I = [I str2double(beta_names{t}(num_idx(1)+4:num_idx(2)-1))];

            num_idx = strfind(beta_names{t},'P_');
            P = [P str2double(beta_names{t}(num_idx+2))];

            num_idx = strfind(beta_names{t},'R_');
            R = [R str2double(beta_names{t}(num_idx+2))];

            num_idx = strfind(beta_names{t},'V');
            V = [V str2double(beta_names{t}(num_idx+1))];
        end
        run_idx = [run_idx; ones(nTrls,1)*run];
        clear SPM beta_names nifti_names
        data = cat(1,data,betas); clear betas %
    end

    X = data; nTrls = length(V); clear data;
    congruency = D==I;

    % Mean-center per run
    for r = 1:nRuns       
         X(run_idx==r,:) = X(run_idx==r,:)-mean(X(run_idx==r,:),1);
         V(run_idx==r) = V(run_idx==r)-mean(V(run_idx==r));
         R(run_idx==r) = R(run_idx==r)-mean(R(run_idx==r)); 
    end

    % Equalize vividness distribution per congruency per run to 



    %% Do the analyses per ROI
    for r = 1:nROIs

        % select ROI
        x = X(:,atlas==r);

        % remove nans
        nan_idx = any(isnan(x));
        x(:,nan_idx) = [];  

        % Get vividness pattern per congruency - cross validate and then
        % downsample test set 
        x_all = []; v_all = []; rj_all = []; rns_all = []; % collect all data for decoding
        for c = 1:2
            idx = congruency==c-1;
            run_idx_cong = run_idx(idx);
            x_cong       = x(idx,:);
            V_cong       = V(idx);
            R_cong       = R(idx);
            RI           = unique(run_idx_cong);

            % Downsample vividness per RJ per run
            x_ds = []; viv = []; rs = []; rns = [];
            for rn = 1:length(RI)
                run = RI(rn);
                rnidx = run_idx_cong==run;

                runx         = x_cong(rnidx,:);
                runv         = V_cong(rnidx);
                runr         = R_cong(rnidx);
                runrn        = run_idx_cong(rnidx);
                [~,~,labels] = unique(runv);
                rj_lab       = unique(runr);

                % log distribution per RJ
                if length(rj_lab)==2 
                v_vals = unique(labels);
                v_distr = nan(2,length(v_vals));
                v_distr_idx = cell(2,length(v_vals));
                for rj = 1:2
                    rjrnidx = find(runr==rj_lab(rj))';
                    rj_labels = labels(rjrnidx);
                    for v = 1:length(v_vals)
                        v_distr_idx{rj,v} = rjrnidx(rj_labels==v_vals(v));
                        v_distr(rj,v) = sum(rj_labels==v_vals(v));
                    end
                end

                % remove vividness ratings with 0 observations
                if any(v_distr(:)==0)
                    if any(v_distr(1,:)==0)
                        v0_idx = find(v_distr(1,:)==0);
                        if any(v_distr(2,:)==0) % check other dim too
                            v0_idx = [v0_idx find(v_distr(2,:)==0)];
                        end
                    elseif any(v_distr(2,:)==0)
                        v0_idx = find(v_distr(2,:)==0);
                        if any(v_distr(1,:)==0) % check other dim too
                            v0_idx = [v0_idx find(v_distr(1,:)==0)];
                        end
                    end
                    v_distr(:,v0_idx) = []; % delete
                    v_distr_idx(:,v0_idx) = [];
                    v_vals = 1:size(v_distr,2);
                end

                % equate distribution by randomly removing trials from
                % largest class 
                new_v_distr_idx = cell(2,length(v_vals));
                for v = 1:length(v_vals)
                    if v_distr(1,v)>v_distr(2,v) % remove from rj=0
                        nVtrls = v_distr(2,v); % smallest trial number
                        larger_idx = v_distr_idx{1,v}(randperm(v_distr(1,v))); % shuffle idx
                        new_v_distr_idx{1,v} = larger_idx(1:nVtrls); % select equal number
                        new_v_distr_idx{2,v} = v_distr_idx{2,v}; % remains the same
                    elseif v_distr(2,v)>v_distr(1,v) % remove from rj=1
                        nVtrls = v_distr(1,v); % smallest trial number
                        larger_idx = v_distr_idx{2,v}(randperm(v_distr(2,v))); 
                        new_v_distr_idx{2,v} = larger_idx(1:nVtrls); 
                        new_v_distr_idx{1,v} = v_distr_idx{1,v};
                    else
                        new_v_distr_idx(:,v) = v_distr_idx(:,v);
                    end
                end
                
                % select downsampled trials
                x_ds = cat(1,x_ds,runx(cell2mat(new_v_distr_idx(:)),:));
                viv = cat(1,viv,runv(cell2mat(new_v_distr_idx(:)))');
                rs = cat(1,rs,runr(cell2mat(new_v_distr_idx(:)))');
                rns = cat(1,rns,runrn(cell2mat(new_v_distr_idx(:))));
                end
            end

            % collect all data
            x_all = cat(1,x_all,x_ds);
            v_all = cat(1,v_all,viv);
            rj_all = cat(1,rj_all,rs);
            rns_all = cat(1,rns_all,rns);

            % cross-validation
            nTrls = length(rs);
            yhat = nan(nTrls,1);
            for rn = 1:length(RI)
                run = RI(rn);
                
                test_idx = find(rns==run);
                train_idx = setdiff(1:nTrls,test_idx);  

                % train             
                beta_ridge = (x_ds(train_idx,:)' * x_ds(train_idx,:) + lambda * eye(size(x_ds(train_idx,:), 2))) \ (x_ds(train_idx,:)' * viv(train_idx));  % Ridge regression coefficients
               
                % get predicted
                yhat(test_idx) = x_ds(test_idx,:)*beta_ridge;
            end

            % calculate accuracy
            acc(sub,r,c,1) = corr(yhat,viv,'Type','Spearman');
            acc(sub,r,c,2) = corr(yhat,rs,'Type','Spearman');          
        end

        % decode vividness and decode RJ over all trials
        nTrls = length(rj_all);
        yhat_V = nan(nTrls,1); yhat_RJ = nan(nTrls,1);
        RI     = unique(rns_all);
        for rn = 1:length(RI)
            run = RI(rn); 

            test_idx = find(rns_all==run);
            train_idx = setdiff(1:nTrls,test_idx);

            % train
            beta_viv = (x_all(train_idx,:)' * x_all(train_idx,:) + lambda * eye(size(x_all(train_idx,:), 2))) \ (x_all(train_idx,:)' * v_all(train_idx));  % Ridge regression coefficients
            beta_rj = (x_all(train_idx,:)' * x_all(train_idx,:) + lambda * eye(size(x_all(train_idx,:), 2))) \ (x_all(train_idx,:)' * rj_all(train_idx));

            % get predicted
            yhat_V(test_idx) = x_all(test_idx,:)*beta_viv;
            yhat_RJ(test_idx) = x_all(test_idx,:)*beta_rj;
        end

        % calculate MSE
        MSE(sub,r,1) = corr(yhat_V,v_all,'Type','Spearman'); %mean((v_all-yhat_V).^2);
        MSE(sub,r,2) = corr(yhat_RJ,rj_all,'Type','Spearman'); %mean((rj_all-yhat_RJ).^2);
    end   
    clear X; 
end

%% Plot the results
figure;
kinds = {'vividness','reality judgement'};
for a = 1:2
    subplot(2,1,a)
    M = squeeze(mean(acc(:,:,:,a)))/sqrt(nSub);
    SEM = squeeze(std(acc(:,:,:,a)));
    barwitherr(SEM,M);
    set(gca,'XTickLabels',cfg.ROIs);
    title(kinds{a}); ylabel('r');
end

figure;
barwitherr(squeeze(std(MSE))./sqrt(nSub),squeeze(mean(MSE)));
set(gca,'XTickLabels',cfg.ROIs);
legend('vividness','reality judgement');
title('Total mean squared error MSE')