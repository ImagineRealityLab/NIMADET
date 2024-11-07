function ModelFittingBehaviour(cfg)

subjects = cfg.subjects(which_subjects);
nSubs    = length(subjects);
nRuns    = 4;

AUC0     = nan(nSubs,2); 
AUC1     = nan(nSubs,2);
BIC0     = nan(nSubs,1); 
BIC1     = nan(nSubs,1);
params   = nan(nSubs,2,2);

P = cell(nSubs,3); % observed, m0, m1
V = cell(nSubs,3);

theta         = [-1 -2]; % starting value P and I

for sub = 1:nSubs
    
    %% Get the behavioural data
    cfg.behDir = fullfile(fileparts(fileparts(cfg.root)),'Results',subjects{sub},'Regressors','Behaviour_matrix');   
    cfg.nRuns  = nRuns;
    [RJ,Vt,Vm,Dp,~,Cond,Pres] = getBehaviour(cfg);
    P{sub,1} = RJ; V{sub,1} = Vt;

    %% Do the model fitting
    fprintf('Model fitting subject %d \n',sub)

    % fit H0
    [params(sub,1,:),~] = fitPRMmodel(RJ,Vt,Vm,Cond,Pres,Dp,theta,1);
    [BIC0(sub), AUC0(sub,:), P{sub,2}, V{sub,2}] = modelPredictions(RJ,Vt,Vm,Cond,Pres,Dp,params(sub,1,:),1,0); 

    % fit H1
    [params(sub,2,:),~] = fitPRMmodel(RJ,Vt,Vm,Cond,Pres,Dp,theta,2);
    [BIC1(sub), AUC1(sub,:), P{sub,3}, V{sub,3}] = modelPredictions(RJ,Vt,Vm,Cond,Pres,Dp,params(sub,2,:),2,0);   

end

%% Save the fit and parameters
outputDir = fullfile(cfg.root,'Results','GroupResults',cfg.outDir);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

save(fullfile(outputDir,'BehaviouralModelFit'));

%% Plot goodness of fit
if cfg.plotGOF
    figure;
    subplot(1,2,1);
    AUC = cat(2,mean(AUC0,2),mean(AUC1,2)); % average AUC over vividness and RJ's
    barwitherr(squeeze(std(AUC))'./sqrt(nSubs),squeeze(mean(AUC))');
    set(gca,'XTickLabels',{'H0 model','H1 model'}); ylim([0.5 0.8])
    ylabel('AUC');
    hold on; scatter((1:2)+(randn(nSubs,1)./10),AUC,60,'filled','k','MarkerFaceAlpha',0.5)

    subplot(1,2,2);
    BIC = [BIC0 BIC1];
    barwitherr(std(BIC)./sqrt(nSubs),mean(BIC));
    set(gca,'XTickLabels',{'H0 model','H1 model'}); %ylim([225 285]); ylabel('BIC')
    fprintf('Mean BIC difference: %.2f \n',mean(BIC0-BIC1))
    ylabel('BIC')
    hold on; scatter((1:2)+(randn(nSubs,1)./10),BIC,60,'filled','k','MarkerFaceAlpha',0.5)
end

% guidelines say between 6 and 10 difference, the model is strongly
% preferred

%% Plot model fit on data
if cfg.plotMF
PP       = nan(nSubs,2,2); % perceptual presence
PPm0     = nan(nSubs,2,2);
PPm1     = nan(nSubs,2,2);

IV       = nan(nSubs,2,2); % imagery vividness
IVm0     = nan(nSubs,2,2);
IVm1     = nan(nSubs,2,2);

for sub = 1:nSubs

    % Get the behavioural data
    cfg     = [];
    cfg.behDir = fullfile(fileparts(fileparts(root)),'Results',subjects{sub},'Regressors','Behaviour_matrix');
    cfg.nRuns  = nRuns;
    [RJ,Vt,Vm,Dp,Cr,Cond,Pres] = getBehaviour(cfg); 

    for c = 1:2
        for p = 1:2
            idx = Cond==c & Pres==(p-1);

            PP(sub,c,p)   = mean(P{sub,1}(idx));
            PPm0(sub,c,p) = mean(P{sub,2}(idx));
            PPm1(sub,c,p) = mean(P{sub,3}(idx));            

            tmp = V{sub,1};%zscore(V{sub,1});
            IV(sub,c,p)   = mean(tmp(idx));

            tmp = V{sub,2};%zscore(V{sub,2});
            IVm0(sub,c,p) = mean(tmp(idx));

            tmp = V{sub,3};%zscore(V{sub,3});
            IVm1(sub,c,p) = mean(tmp(idx));
        end
    end

end

% plot condition effects
names = {'PP','PPm0','PPm1','IV','IVm0','IVm1'};
figure;
for i = 1:6
    data = eval(names{i});
    M    = squeeze(mean(data));
    SEM  = squeeze(std(data))./sqrt(nSubs);

    subplot(2,3,i)
    barwitherr(SEM',M');
    title(names{i})
    if contains(names{i},'PP'); ylim([0 1]); else; ylim([1 4]); end

    set(gca,'XTickLabels',{'Stim abs','Stim pres'})
end

% plot model fit on top of oberved data
names = {'PP','PPm0','PPm1';'IV','IVm0','IVm1'};
figure; cs = {'y','c'};
for i = 1:2
    data = eval(names{i,1});
    M    = squeeze(mean(data));
    SEM  = squeeze(std(data))./sqrt(nSubs);

    subplot(1,2,i)
    b = barwitherr(SEM',M');
    title(names{i})

    set(gca,'XTickLabels',{'Stim abs','Stim pres'})

    % add model fit
    for m = 1:2
        data = eval(names{i,1+m});
        M    = squeeze(mean(data));
        SEM  = squeeze(std(data))./sqrt(nSubs);
        for j = 1:2
            hold on; errorbar(b(j).XData+b(j).XOffset,M(j,:),SEM(j,:),'o','Color',cs{m},'LineWidth',2);
        end

    end
end
end
