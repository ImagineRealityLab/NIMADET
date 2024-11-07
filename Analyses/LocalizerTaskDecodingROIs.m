function LocalizerTaskDecodingROIs2(cfg)

nVox     = 200;
nPerm    = 25;

nSub = length(cfg.subjects);


%cfg.ROI   = fullfile('D:\NIMADET\Analyses\ROI_masks\Conjunction\conjunction_atlasFWE.nii');
[~,atlas] = read_nii(cfg.ROI);
mask      = atlas>0; atlas = atlas(mask);
nROIs     = length(unique(atlas(:)));

trial_type{1} = 'congruency == 1 & P == 1';
trial_type{2} = 'congruency == 1 & P == 0';
trial_type{3} = 'congruency == 0 & P == 1';
trial_type{4} = 'congruency == 0 & P == 0';
nTT = length(trial_type);

Acc   = nan(nSub,nROIs,nTT); 
Dat   = []; % save classifier evidence per trial

if cfg.permutation
    PAcc = nan(nSub,nROIs,nTT,nPerm);
end

rng(1,'twister'); % for reproducability

for sub = 1:nSub
    fprintf('Processing subject %s \n',  cfg.subjects{sub})

    % out dir
    outDir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.outDir);
    if ~exist(outDir,'dir'); mkdir(outDir); end

    if ~exist(fullfile(outDir,'ROI_accuracy.mat'),'file')

        %% prepare the data
        % get the t-map and mask
        [~,tmap] = read_nii(fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.stim));
        tmap     = tmap(mask(:));
       
        % task data
        dataDir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.data{1});
        nRuns   = length(dir(dataDir))-2;
        D = []; I = []; P = []; R = []; V = []; data = []; run_idx = [];
        for run = 1:nRuns

            fprintf('Reading in task betas for run %d... \n',run)

            beta_niftis = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),'trl_*.nii');
            nifti_names = cell(length(beta_niftis),1); % re-sort
            beta_names  = cell(length(nifti_names),1);
            for n = 1:length(beta_niftis)
                nifti_names{n} = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),sprintf('trl_%d_*.nii',n));
                [~,beta_names{n}] = fileparts(nifti_names{n});
            end

            nTrls = length(beta_names);

            % remove blocks failed imagery check
            if cfg.ima_check
                load(fullfile(cfg.root,'Data',cfg.subjects{sub},...
                    'Behaviour',sprintf('MT_%s_run%d.mat',cfg.subjects{sub},run)),'C')
                if any(C==0)
                    task_ridx = [];
                    for c = 1:length(C)
                        idx = (c-1)*(nTrls/length(C))+1:c*(nTrls/length(C));
                        if C(c)==0
                            task_ridx = [task_ridx,idx];
                        end
                    end
                    beta_names(task_ridx) = [];
                    nifti_names(task_ridx) = [];
                end
            end

            nTrls = length(beta_names);
            betas = nan(nTrls,sum(mask(:)));
            for t = 1:nTrls
                [~,Y] = read_nii(nifti_names{t});
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
            data = cat(1,data,betas); clear betas
        end

        % define congruency
        congruency = D == I;
        task_oris = unique(D);

        % loop over trial types
        task_data = cell(nTT,1); Y_task = cell(nTT,1); 
        task_ridx = cell(nTT,1); R_task = cell(nTT,1); % responses
        for tt = 1:nTT

            trial_type_idx = eval(trial_type{tt});

            task_data{tt} = data(trial_type_idx,:);
            Y_task{tt}    = D(trial_type_idx);
            task_ridx{tt} = run_idx(trial_type_idx);
            R_task{tt}    = [V(trial_type_idx)' R(trial_type_idx)'];

            % first balance, then z-score per run
            X = []; Y = []; Res = []; RN = []; 
            RI = unique(task_ridx{tt}); nRuns = length(RI);
            for r= 1:nRuns

                run = RI(r);

                xrun = task_data{tt}(task_ridx{tt}==run,:);
                yrun = double(Y_task{tt}(task_ridx{tt}==run)==task_oris(1));
                rrun = R_task{tt}(task_ridx{tt}==run,:);

                % balance
                idx = balance_trials(yrun'+1,'downsample');
                xrun = xrun(cell2mat(idx),:);
                yrun = yrun(cell2mat(idx))';
                rrun = rrun(cell2mat(idx),:);

                Y  = cat(1,Y,yrun);
                RN = cat(1,RN,ones(length(yrun),1)*run);
                Res  = cat(1,Res,rrun);

                % mean centre per run
                X = cat(1,X,xrun-nanmean(xrun,2)); % mean over voxels
            end
            Y_task{tt} = Y;
            task_ridx{tt} = RN;
            task_data{tt} = X;
            R_task{tt} = Res;
        end

        % localizer data
        dataDir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.data{2});
        nRuns   = length(dir(dataDir))-2; run_idx = [];
        X       = cell(nRuns,1); Y = cell(nRuns,1);
        for run = 1:2
            run_dir = fullfile(dataDir,sprintf('RUN_%02d',run));
            if contains(dataDir,'Poldrack')
                trl_niftis = str2fullfile(run_dir,'*.nii');
                nTrls = length(trl_niftis); X{run} = nan(nTrls,sum(mask(:))); Y{run} = nan(nTrls,1);
                for t = 1:nTrls
                    [~,beta] = read_nii(trl_niftis{t});
                    X{run}(t,:) = beta(mask);
                    [~,name,~] = fileparts(trl_niftis{t});
                    num_idx = strfind(name,'stm_');
                    Y{run}(t) = str2double(name(num_idx+4:end));
                end

            else
                load(fullfile(run_dir,'SPM.mat'),'SPM');
                bidx = find(contains(SPM.xX.name,'bf(1)') & contains(SPM.xX.name,'trl'));

                nTrls = length(bidx); X{run} = nan(nTrls,sum(mask(:))); Y{run} = nan(nTrls,1);
                for t = 1:nTrls
                    [~,beta] = read_nii(fullfile(run_dir,...
                        sprintf('beta_%04d.nii',bidx(t))));
                    X{run}(t,:) = beta(mask);
                    num_idx = [strfind(SPM.xX.name{bidx(t)},'stm_') strfind(SPM.xX.name{bidx(t)},'*bf(1)')];
                    Y{run}(t)   = str2double(SPM.xX.name{bidx(t)}(num_idx(1)+4:num_idx(2)-1));
                    clear beta
                end
            end
            run_idx = [run_idx;ones(nTrls,1)*run];
        end

       % only select task oris
        if numel(X)==2
            X = cat(1,X{1},X{2}); Y = cat(1,Y{1},Y{2});
        else
            X = X{1}; Y = Y{1};
        end
        run_idx = run_idx(ismember(Y,task_oris),:);
        X = X(ismember(Y,task_oris),:);
        Y = Y(ismember(Y,task_oris),:);
        
        % mean-centre per run
        for run = 1:nRuns
            X(run_idx==run,:) = X(run_idx==run,:)-nanmean(X(run_idx==run,:),2); % mean over voxels
        end

        loc_data = X; loc_labels = Y; loc_ridx = run_idx; 
        Y_loc = loc_labels; Y_loc = double(loc_labels==task_oris(1));
        clear X Y run_idx

        %% Do the decoding per ROI
        acc      = nan(nROIs,nTT); % train loc - train task
        pAcc     = nan(nROIs,nTT,nPerm);
        dat      = [];

        for r = 1:nROIs

            fprintf('Decoding ROI %d out of %d \n',r,nROIs)

            for i = 1:nTT

                x{1} = loc_data(:,atlas==r); x{2} = task_data{i}(:,atlas==r);

                % get n most responsive voxels
                if size(x{1},2) > nVox
                    [~,act_ind] = sort(tmap(atlas==r),'descend');
                    act_ind = act_ind(1:nVox);
                    x{1} = x{1}(:,act_ind);
                    x{2} = x{2}(:,act_ind);
                end

                % remove nans
                nan_idx = any(isnan(x{1})) | any(isnan(x{2}));
                x{1}(:,nan_idx) = []; x{2}(:,nan_idx) = [];

%                 % mean center per run
%                 RN = task_ridx{i};
%                 for rn = 1:length(unique(RN))
%                    x{2}(RN==rn,:) = x{2}(RN==rn,:)-mean(x{2}(RN==rn,:),1);
%                 end

                % MEAN CENTER OVER VOXELS
                x{1} = x{1}-mean(x{1},2);
                x{2} = x{2}-mean(x{2},2);

                % classify
                if ~isempty(x{1})

                    if length(unique(Y_task{i}))>1 % need to have at least 1 trial per orientation

                        % LDA decoding
                        cfgD = []; cfgD.gamma = cfg.gamma; 
                        cfgD.demean = 'no';
                        decoder = train_LDA(cfgD,Y_loc,x{1}');
                        yhat = decode_LDA(cfgD,decoder,x{2}');
                        acc(r,i) = mean((yhat>0)==Y_task{i}');

                        % save classifier evidence
                        evi = yhat; evi(Y_task{i}==0)=yhat(Y_task{i}==0)*-1;
                        nTrls = length(yhat);

                        tmp = [ones(nTrls,1).*sub ones(nTrls,1).*r ones(nTrls,1).*i Y_task{i} R_task{i} yhat' evi'];
                        dat = cat(1,dat,tmp);

                        if cfg.permutation
                            for p = 1:nPerm                                
                                permY = Y_loc(randperm(length(Y_loc))); 
                                cfgD = []; cfgD.gamma = cfg.gamma; 
                                cfgD.demean = 'no';
                                decoder = train_LDA(cfgD,permY,x{1}');
                                yhat = decode_LDA(cfgD,decoder,x{2}');
                                pAcc(r,i,p) = mean((yhat>0)==Y_task{i}');
                            end
                        end
                       
                    else
                        fprintf('\t Only 1 orientation for idx %d congruency %d \n',t,c)

                    end

                end
            end
        end

        % write the results
        if cfg.permutation
            save(fullfile(outDir,'ROI_accuracy'),'acc','pAcc','dat');
        else
            load(fullfile(outDir,'ROI_accuracy'),'acc','dat');
        end

    else
        fprintf('\t Decoding already done \n')
        load(fullfile(outDir,'ROI_accuracy'),'acc','pAcc','dat');
        if ~exist('dat','var'); dat = []; end
    end

    % concatenate subs
    Acc(sub,:,:) = acc;
    Dat   = cat(1,Dat,dat);

    if cfg.permutation
        PAcc(sub,:,:,:) = pAcc(:,:,:);
    end

end

%% Save group results
outDir = fullfile(cfg.root,'Results','GroupResults',cfg.outDir);
if ~exist(outDir,'dir'); mkdir(outDir); end

if ~exist(fullfile(outDir,'ROI_accuracy'),'file')
    save(fullfile(outDir,'ROI_accuracy'),'Dat','Acc','PAcc');
else
    load(fullfile(outDir,'ROI_accuracy'),'Dat','Acc','PAcc');
end

%% Stats
if cfg.permutation

    pVals = nan(nROIs,nTT);
    nBtstrp = 10000;

    for r = 1:nROIs
        for tt = 1:nTT
            bAcc1 = nan(nBtstrp,1);
            for b = 1:nBtstrp
                tmp1 = nan(nSub,1);
                for sub = 1:nSub
                    tmp1(sub) = PAcc(sub,r,tt,randi(nPerm));
                end
                bAcc1(b) = mean(tmp1);
            end

            pVals(r,tt) = sum(bAcc1>squeeze(mean(Acc(:,r,tt))))./nBtstrp;
        end
    end
    save(fullfile(outDir,'ROI_accuracy'),'pVals','-append');
end

%% Plotting
ROI_names = {'ROI'};
for r = 1:nROIs
    figure(1); % decoding accuracy
    subplot(2,ceil(nROIs/2),r);
    dat = squeeze(Acc(:,r,:));
    M   = mean(dat);
    SEM = std(dat)./sqrt(nSub);
    b = barwitherr(SEM,M);
    hold on; plot(xlim,[0.5 0.5],'k--');
    ylim([0.45 0.55]);
    %hold on; scatter(b.XData+(randn(nTT,nSub)./10)',dat,40,'k','filled','MarkerFaceAlpha',0.5)
    set(gca,'XTickLabels',{'cong pres','cong abs','inco pres','inco abs'});
    title(ROI_names{r});
end

%% Look at influence responses on decodeability
% subject     = Dat(:,1);
% roi         = Dat(:,2);
% condition   = Dat(:,3);
% orientation = Dat(:,4);
% vividness   = Dat(:,5);
% response    = Dat(:,6);
% distance    = Dat(:,7);
% accuracy    = Dat(:,8);
% 
% %tmp = Dat;
% %tmp(:,1) = Dat(:,1); tmp(:,2) = ones(length(tmp),1); tmp(:,3:8) = Dat(:,2:7); Dat = tmp; clear tmp
% idx  = roi==1 & condition==2;
% data = table(subject(idx),condition(idx),orientation(idx),vividness(idx),response(idx),distance(idx),accuracy(idx));
% data.Properties.VariableNames  = {'subject','condition','orientation','vividness','response','distance','accuracy'};
% lme = fitlme(data,'accuracy~vividness+response+(1|vividness)+(1|response)+(1|subject)');
% lme = fitlme(data,'accuracy~vividness+(1|vividness)+(1|subject)');

beh_betas = nan(nSub,nROIs,2,5);
for sub = 1:nSub

    for r = 1:nROIs

        dat = squeeze(Dat(Dat(:,1)==sub & Dat(:,2)==r,3:end));
        X   = [ones(length(dat),1) dat(:,3) dat(:,4)];
        b   = X\dat(:,6);
        beh_betas(sub,r,:,1) = b(2:end);

        for tt = 1:nTT
            d = squeeze(dat(dat(:,1)==tt,:));
            X   = [ones(length(d),1) d(:,3) d(:,4)];
            b   = X\d(:,6);
            beh_betas(sub,r,:,tt+1) = b(2:end);
        end
    end
end


M = nan(nROIs,2,5); SEM = nan(nROIs,2,5);
for r = 1:nROIs
    for re = 1:2
        for tt = 1:nTT+1
            dat = beh_betas(:,r,re,tt);
            dat(isnan(dat)) = [];
            M(r,re,tt) = mean(dat); SEM(r,re,tt) = std(dat)./sqrt(length(dat));
        end
    end
end

figure;
for r = 1:nROIs
    subplot(2,ceil(nROIs/2),r);
    b = barwitherr(squeeze(SEM(r,:,:)),squeeze(M(r,:,:))); %legend('Overall','Cong pres','Cong abs','Inco pres','Inco abs');
%     for i = 1:nTT+1
%         dat = squeeze(beh_betas(:,r,:,i));
%         hold on; scatter(b(i).XData+b(i).XOffset+(randn(2,nSub)./50)',dat,40,'k','filled','MarkerFaceAlpha',0.5)
%     end
    set(gca,'XTicklabels',{'Vividness','Detection'});
    ylabel('Betas')
    title(ROI_names{r})
end
nan_idx = any(any(any(beh_betas==0,4),3),2);
[h,p]   = ttest(squeeze(beh_betas(~nan_idx,:,:,:)));

%% Split pres/abs and low/high vividness
acc_R = nan(nSub,nROIs,5,2,2);
for sub = 1:nSub
    for r = 1:nROIs
        dat = squeeze(Dat(Dat(:,1)==sub & Dat(:,2)==r,3:end));
        R   = [dat(:,3)>(median(dat(:,3))) dat(:,4)];
        for i1 = 1:2 % vividness vs response
            for i2 = 1:2 % low versus high
                idx = R(:,i1)==(i2-1);
                acc_R(sub,r,1,i1,i2) = mean((dat(idx,5)>0)==dat(idx,2));
            end
        end

        % per trial type
        for tt = 1:nTT
            d = dat(dat(:,1)==tt,2:end);
            re = [d(:,2)>(median(d(:,2))) d(:,3)];
            for i1 = 1:2 % vividness vs response
                for i2 = 1:2 % low vs high
                    idx = re(:,i1)==(i2-1);
                    acc_R(sub,r,tt+1,i1,i2) = mean((d(idx,4)>0)==d(idx,1));
                end
            end
        end
    end
end

M = nan(nROIs,5,2,2); SEM = nan(nROIs,5,2,2);
for r= 1:nROIs
    for tt = 1:5
        for i1 = 1:2
            for i2 = 1:2
                tmp = squeeze(acc_R(:,r,tt,i1,i2));
                tmp(isnan(tmp)) = [];
                M(r,tt,i1,i2) = mean(tmp); SEM(r,tt,i1,i2) = std(tmp)./sqrt(length(tmp));
            end
        end
    end
end

figure; names = {'all','cong pres','cong abs','inco pres','inco abs'};
if nROIs > 1
for r = 1:nROIs
    for tt = 1:5
        subplot(5,12,r+(tt-1)*12)
        barwitherr(squeeze(SEM(r,tt,:,:)),squeeze(M(r,tt,:,:))); hold on; plot(xlim,[0.5 0.5],'k--');
        ylim([0.4 0.6]); set(gca,'XTickLabels',{'Vividness','Response'});
        %legend('Low','High'); 
        title(names{tt});
    end
end
else
    r = 1;
    for tt = 1:5
        subplot(1,5,tt)
        barwitherr(squeeze(SEM(r,tt,:,:)),squeeze(M(r,tt,:,:))); hold on; plot(xlim,[0.5 0.5],'k--');
        ylim([0.4 0.6]); set(gca,'XTickLabels',{'Vividness','Response'});
        legend('Low','High'); 
        title(names{tt});
    end
end
%% Check individual differences
cfg.plotIndividual = false; 
cfg.plotting = false;

[Cr,Dp,FA,H,V] = BehaviourAnalysis(cfg); % per subject now 
shift = Cr(:,2)-Cr(:,1);
Dp = mean(Dp,2); Cr = mean(Cr,2);
V = squeeze(mean(mean(V,3),2));

vals = {'Dp','Cr','shift','V'};
if nROIs == 1
figure;
for v = 1:length(vals)
    subplot(2,2,v)
    m = eval(vals{v});
    dat1 = squeeze(Acc(m<median(m),:,:));
    dat2 = squeeze(Acc(m>median(m),:,:));
    SEM  = [std(dat1)./sqrt(sum(m<median(m))); std(dat2)./sqrt(sum(m>median(m)))];
    M    = [mean(dat1); mean(dat2)];
    barwitherr(SEM',M'); ylim([0.45 0.6]);
    hold on; plot(xlim,[0.5 0.5],'k--')
    title(vals{v})
    set(gca,'XTickLabels',{'CP','CA','IP','IA'})
end

else
    for v= 1:length(vals)

        m = eval(vals{v});
        figure(v);
        for r = 1:nROIs 
            subplot(2,6,r)
            dat1 = squeeze(Acc(m<median(m),r,:));
            dat2 = squeeze(Acc(m>median(m),r,:));
            SEM  = [std(dat1)./sqrt(sum(m<median(m))); std(dat2)./sqrt(sum(m>median(m)))];
            M    = [mean(dat1); mean(dat2)];
            barwitherr(SEM',M'); ylim([0.45 0.6]);
            hold on; plot(xlim,[0.5 0.5],'k--')
            title(ROI_names{r})
            set(gca,'XTickLabels',{'CP','CA','IP','IA'})
        end
    end
end
%% Calculate accuracy for different viv levels
vivAcc = nan(nSub,4); 
for sub = 1:nSubs
    for v = 1:4
        idx = Dat(:,1)==sub & Dat(:,5)==v;
        if sum(idx)>20
            vivAcc(sub,v) = mean((Dat(idx,7)>0)==Dat(idx,4));
        end
    end
end

M = nan(1,4); SEM = nan(1,4);
for v = 1:4
    dat = vivAcc(:,v);
    dat(isnan(dat)) = [];
    M(v) = mean(dat); SEM(v) = std(dat)./sqrt(length(dat));
end

figure;
barwitherr(SEM,M); ylim([0.45 0.55]); 