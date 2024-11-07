function ExtractModelRS(cfg)

subjects = cfg.subjects(which_subjects);
nSubs    = length(subjects);
nRuns    = 4;

outputDir = fullfile(cfg.root,'Results','GroupResults',cfg.outDir);
load(fullfile(outputDir,'BehaviouralModelFit'),'params');

%% Simulate RS for neuroimaging analysis
RS_percond = nan(nSubs,2,2); % congruency x presence
RS_perresp = nan(nSubs,4); % per response type
for sub = 1:nSubs
    
    %% Get the behavioural data
    cfg.behDir = fullfile(fileparts(fileparts(cfg.root)),'Results',subjects{sub},'Regressors','Behaviour_matrix');   
    cfg.nRuns  = nRuns;
    [RJ,Vt,Vm,Dp,~,Cond,Pres] = getBehaviour(cfg);

    %% Simulate signal with fitted parameters
    [~, ~, RJ_pred, v_pred,RS_sim] = modelPredictions(RJ,Vt,Vm,Cond,Pres,Dp,params(sub,2,:),2,0); 
    RJ_pred = double(RJ_pred>0.5); RS_sim = zscore(RS_sim); 

    %% Get RS per condition
    for c = 1:2
        for pres = 1:2
            idx = Cond==c & Pres==(pres-1);
            RS_percond(sub,c,pres) = mean(RS_sim(idx));
        end        
    end

    % per response type
    count = 1;
    for rj = 1:2
        idx = RJ_pred==(rj-1);
        tmp = RS_sim(idx);
        V_split = double(v_pred(idx)>median(v_pred(idx)));
        for v = 1:2
            RS_perresp(sub,count) = mean(tmp(V_split==(v-1)));
            count = count+1;
        end
    end

end

% plot the model RS
figure;
subplot(2,1,1);
M = squeeze(mean(RS_percond)); SEM = squeeze(std(RS_percond))./sqrt(nSubs);
b = barwitherr(SEM',M');
for c = 1:2
    hold on; s = scatter((b(c).XData+b(c).XOffset)+randn(nSubs,1)./40,...
        squeeze(RS_percond(:,c,:)),30,'filled','k','MarkerFaceAlpha',0.5);
end
set(gca,'XTickLabel',{'Absent','Present'}); legend('Congruent','Incongruent')

subplot(2,1,2);
nan_idx = squeeze(any(isnan(RS_perresp),2));
RS_perresp(nan_idx,:) = [];
M = squeeze(mean(RS_perresp)); SEM = squeeze(std(RS_perresp))./sqrt(size(RS_perresp,1));
barwitherr(SEM,M)
set(gca,'XTickLabel',{'Congruent','Incongruent'});
