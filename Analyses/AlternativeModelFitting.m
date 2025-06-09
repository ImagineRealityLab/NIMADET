%ddpath(genpath('D:\NIMADET\ModelSimulations\Analyses'))
root = 'D:\NIMADET\ModelSimulations\NewSimulations\Revision1_Neuron';
cd(root)
addpath('Models');

subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11',...
    'S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24',...
    'S25','S26','S27','S28','S29','S30','S31','S32','S33','S34','S35','S36'};
which_subjects = [1:3 5:6 9 11:16 18 20 21 23 26:28 31:36]; % 7 has too little vividness variation

subjects = subjects(which_subjects);
nSubs    = length(subjects);
nRuns    = 4;

AUC     = nan(nSubs,4,2); 
BIC     = nan(nSubs,4); 
params   = nan(nSubs,4,3);

P = cell(nSubs,3); % observed, m0, m1
V = cell(nSubs,3);

load('AllModelParameters')
params(6,:,:) = [];

for sub = 1:nSubs
    
    %% Get the behavioural data
    cfg     = [];
    cfg.behDir = fullfile(fileparts(fileparts(fileparts(root))),'Results',subjects{sub},'Regressors','Behaviour_matrix');   
    cfg.nRuns  = nRuns;
    [RJ,Vt,Vm,Dp,Cr,Cond,Pres] = getBehaviour(cfg);
    P{sub,1} = RJ; V{sub,1} = Vt;

    %% Do the model fitting
    fprintf('Model fitting subject %d \n',sub)

    % Null model
    %theta = [-1 -2];
    %[params(sub,1,1:2),~] = fitAlternativeModels(RJ,Vt,Vm,Cond,Pres,Dp,theta,1);
    [BIC(sub,1), AUC(sub,1,:), P{sub,2}, V{sub,2}] = modelPredictions(RJ,Vt,Vm,Cond,Pres,Dp,params(sub,1,:),1,0); 

    % Reality threshold model
    %theta = [-1 -2 1];
    %[params(sub,2,:),~] = fitAlternativeModels(RJ,Vt,Vm,Cond,Pres,Dp,theta,2);
    [BIC(sub,2), AUC(sub,2,:), P{sub,3}, V{sub,3}] = modelPredictions(RJ,Vt,Vm,Cond,Pres,Dp,params(sub,2,:),2,0);  
 
    % Multiplicative model
    %theta = [-1 -2];
    %[params(sub,3,1:2),~] = fitAlternativeModels(RJ,Vt,Vm,Cond,Pres,Dp,theta,3);
    [BIC(sub,3), AUC(sub,3,:), P{sub,4}, V{sub,4}] = modelPredictions(RJ,Vt,Vm,Cond,Pres,Dp,params(sub,3,:),3,0); 

    % Sublinear model
    %theta = [-1 -2 1];
    %[params(sub,4,:),~] = fitAlternativeModels(RJ,Vt,Vm,Cond,Pres,Dp,theta,4);
    [BIC(sub,4), AUC(sub,4,:), P{sub,5}, V{sub,5}] = modelPredictions(RJ,Vt,Vm,Cond,Pres,Dp,params(sub,4,:),4,0);  

end


%% Save the fit and parameters
%save('ModelFitMultiplicativeSublinear','params');

%% Plot goodness of fit
% remove nans
nan_idx = squeeze(any(any(isnan(AUC),2),3));
AUC(nan_idx,:,:) = []; nSubs = size(AUC,1);
BIC(6,:) = [];

figure;
subplot(1,2,1); AUCm = squeeze(mean(AUC,3));
barwitherr(squeeze(std(AUCm))'./sqrt(nSubs),squeeze(mean(AUCm))');
set(gca,'XTickLabels',{'Null','Additive','Multiplicative','Sublinear'}); ylim([0.4 0.8])
ylabel('AUC');
hold on; scatter((1:4)+(randn(nSubs,1)./10),AUCm,60,'filled','k','MarkerFaceAlpha',0.5)

subplot(1,2,2);
BICgroup = sum(BIC);
bar(BICgroup-min(BICgroup))
set(gca,'XTickLabels',{'Null','Additive','Multiplicative','Sublinear'}); 
ylabel('Summed BIC'); title('Bayesian information criterion')

% guidelines say between 6 and 10 difference, the model is strongly
% preferred

%% Plot model fit on data 
PP       = nan(nSubs,2,2); % perceptual presence
PPm      = nan(4,nSubs,2,2);

IV       = nan(nSubs,2,2); % imagery vividness
IVm      = nan(4,nSubs,2,2);

for sub = 1:nSubs

    % Get the behavioural data
    cfg     = [];
    cfg.behDir = fullfile(fileparts(fileparts(fileparts(root))),'Results',subjects{sub},'Regressors','Behaviour_matrix');
    cfg.nRuns  = nRuns;
    [RJ,Vt,Vm,Dp,Cr,Cond,Pres] = getBehaviour(cfg); 

    for c = 1:2
        for p = 1:2
            idx = Cond==c & Pres==(p-1);

            PP(sub,c,p)   = mean(P{sub,1}(idx));
            IV(sub,c,p)   = mean(V{sub,1});

            for m = 1:4
                PPm(m,sub,c,p) = mean(P{sub,m+1}(idx));
                IVm(m,sub,c,p) = mean(V{sub,m+1}(idx));
            end       

        end
    end

end

% plot condition effects
names = {'data','Null','Additive','Multiplicative','Sublinear'};
figure;
for m = 1:5
    if m == 1
        data = PP;
    else
        data = squeeze(PPm(m-1,:,:,:));
    end
    M    = squeeze(mean(data));
    SEM  = squeeze(std(data))./sqrt(nSubs);
    subplot(2,5,m) 
    barwitherr(SEM',M');
    title(names{m}); ylim([0 1]) 
    set(gca,'XTickLabels',{'Stim abs','Stim pres'})

    if m == 1
        data = IV;
    else
        data = squeeze(IVm(m-1,:,:,:));
    end
    M    = squeeze(mean(data));
    SEM  = squeeze(std(data))./sqrt(nSubs);
    subplot(2,5,m+5) 
    barwitherr(SEM',M');
    title(names{m}); ylim([1 4]) 
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
