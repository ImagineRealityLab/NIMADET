function [Cr,Dp,FA,H,V,acc,Vpr] = BehaviourAnalysis(cfg)
% function BehaviourAnalysis(cfg)

nSub = length(cfg.subjects);
FA   = nan(nSub,2); H = nan(nSub,2);
Dp   = nan(nSub,2); Cr = nan(nSub,2);
V    = nan(nSub,2,2); % per condition per presence
Vpr  = nan(nSub,2,2,2); 
rmIC = nan(nSub,1);
acc  = nan(nSub,2);

for sub = 1:nSub
    % output dir
    outputDir = fullfile(cfg.root,'Results',cfg.subjects{sub});
    if ~exist(outputDir,'dir'); mkdir(outputDir); end
    
    % get the response mapping
    load(fullfile(cfg.root,'Data',cfg.subjects{sub},'Behaviour',['RMs_' cfg.subjects{sub} '.mat']),'responseMappings')
    
    % get the data
    data = str2fullfile(fullfile(cfg.root,'Data',cfg.subjects{sub},'Behaviour'),...
        'MT_*.mat');
    nRuns = length(data);
    resps = []; trl_info = []; ima_check = []; blck_info = [];
    for r = 1:nRuns
        
        load(data{r},'R','C','trials','blocks', 'orientations','miniblocks');
        p_ori = miniblocks; clear miniblocks;
        i_ori = reshape(repmat(blocks,1,length(p_ori)/length(blocks))',length(p_ori),1); 

        if responseMappings(r) == 2
            R(:,:,1) = 5-R(:,:,1); % reverse 
        end

        trl_info  = cat(1,trl_info,trials);
        resps     = cat(1,resps,R);
        ima_check = cat(1,ima_check,C);
        blck_info = cat(1,blck_info,[p_ori,i_ori]);
        
        clear R C trials blocks
    end
    
    % remove incorrect imagery blocks
    ima_check = reshape(repmat(ima_check,1,2),length(blck_info),1);  
    idx = ima_check==0;
    blck_info(idx,:) = []; resps(idx,:,:) = [];
    trl_info(idx,:) = []; 
    fprintf('Removed %d mini-blocks due to failed imagery check \n',sum(idx));
    rmIC(sub,1) = sum(idx);
    
    % analyze the responses
    clear idx
    idx{1} = blck_info(:,1) == blck_info(:,2); % imagine and detect same stimulus
    idx{2} = blck_info(:,1) ~= blck_info(:,2);
    
    trl_resps = cell(2,1); trl_pres = cell(2,1);
    for c = 1:2
        % get the right trials
        trl_resps{c} = reshape(permute(resps(idx{c},:,:),[2,1,3]),sum(idx{c})*size(resps,2),size(resps,3));
        trl_pres{c}  = reshape(trl_info(idx{c},:)',sum(idx{c})*size(trl_info,2),1);
        
        % calculate accuracy
        results.acc(c) = mean(trl_resps{c}(:,3)==trl_pres{c});
        acc(sub,c) = results.acc(c);
        
        % calculate FA's and H's
        for r = 1:2
            rps       = sum(trl_resps{c}(:,3) == 1 & trl_pres{c} == (r-1));
            presTrls  = sum(trl_pres{c}==(r-1));
            if rps == 0; rps = 0.5; elseif rps/presTrls == 1; rps = rps-0.5; end
            results.FA_Hs(c,r) = rps/presTrls;
            
            for p = 1:2
                results.V(c,r,p) = nanmean(trl_resps{c}(trl_resps{c}(:,3)==(r-1) & trl_pres{c} == (p-1),1));
                Vpr(sub,c,r,p) = results.V(c,r,p);

                V(sub,c,p) = nanmean(trl_resps{c}(trl_pres{c} == (p-1),1));
            end
        end
        
        % get d-prime
        [results.D(c),results.C(c)] = dprime(results.FA_Hs(c,2),results.FA_Hs(c,1));
    end
    
    % analyze per block
    cong = double(blck_info(:,1)~=blck_info(:,2))+1;
    acc_blck = cell(2,1);
    for c = 1:2
        blck_idx = find(cong==c);
        nBlcks = length(blck_idx);
        acc_blck{c} = nan(nBlcks,1);
        for b = 1:nBlcks
            acc_blck{c}(b) = mean(squeeze(resps(blck_idx(b),:,3))...
                == trl_info(blck_idx(b),:));
        end        
    end
    
    % orientations
    results.orientations = orientations; clear orientations
    
    % save them
    save(fullfile(outputDir,'behResults'),'results')
    
    FA(sub,:) = results.FA_Hs(:,1);
    H(sub,:)  = results.FA_Hs(:,2);
    Dp(sub,:) = results.D;
    Cr(sub,:)  = results.C;
    %V(sub,:,:,:) = results.V;

    % plot individual?
    if cfg.plotIndividual
        
        % plot the things
        figure; % detection performance
        subplot(2,2,1:2)
        b = bar(results.FA_Hs'); b(1).FaceColor = [0 0 0.5]; b(2).FaceColor = [0.5 0 0]; ylim([0 1])
        legend('Congruent','Incongruent'); set(gca,'XTickLabels',{'False Alarms','Hits'});
        subplot(2,2,3); b = bar(1,results.D(1)); b.FaceColor = [0 0 0.5];
        hold on; b = bar(2,results.D(2)); b.FaceColor = [0.5 0 0]; title("d'");
        subplot(2,2,4); b = bar(1,results.C(1)); b.FaceColor = [0 0 0.5];
        hold on; b = bar(2,results.C(2)); b.FaceColor = [0.5 0 0]; title("C");
        
        figure; % vividness
        subplot(2,1,1);
        b = bar(squeeze(results.V(:,:,1))'); b(1). FaceColor = [0 0 0.5];
        b(2).FaceColor = [0.5 0 0]; legend('Congruent','Incongruent');
        set(gca,'XTickLabels',{'CR','FA'}); ylabel('Vividness');
        subplot(2,1,2);
        b = bar(squeeze(results.V(:,:,2))'); b(1). FaceColor = [0 0 0.5];
        b(2).FaceColor = [0.5 0 0]; legend('Congruent','Incongruent');
        set(gca,'XTickLabels',{'Miss','Hit'}); ylabel('Vividness');        
    end
    
    clear results
end

% plot group results
if cfg.plotting
    
    names = {'FAs','Hits',"D'",'C'};
    data  = cat(3,FA,H,Dp,Cr);
    
    figure;
    for d = 1:4
        subplot(2,2,d); b = bar(1,squeeze(mean(data(:,1,d),1)));
        b.FaceColor = [0 0 0.5]; 
        hold on; scatter(1+randn(nSub,1)./10,squeeze(data(:,1,d)),40,...
            'filled','MarkerFaceColor','k','MarkerFacealpha',0.5)
        hold on; b = bar(2,squeeze(mean(data(:,2,d),1)));
        b.FaceColor = [0.5 0 0]; 
        hold on; scatter(2+randn(nSub,1)./10,squeeze(data(:,2,d)),40,...
            'filled','MarkerFaceColor','k','MarkerFacealpha',0.5)
        title(names{d})
    end

    % plot influence stim and condition on RJ and vividness 
    figure;
    subplot(2,1,1);
    RJ = cat(3,FA,H);
    M = squeeze(mean(RJ,1)); SEM = squeeze(std(RJ,[],1))./sqrt(nSub);
    b = barwitherr(SEM',M'); 
    b(1).FaceColor = [0 0 0.5];
    b(2).FaceColor = [0.5 0 0]; 
    legend('Congruent','Incongruent'); ylim([0 1])
    set(gca,'XTickLabels',{'Absent','Present'}); ylabel('Detection response');
    for p = 1:2
        hold on; scatter([p-0.15 p+0.15]+randn(nSub,1)./20,squeeze(RJ(:,:,p)),40,...
            'filled','MarkerFaceColor','k','MarkerFacealpha',0.5)
    end

    subplot(2,1,2);
    M = squeeze(mean(V,1)); SEM = squeeze(std(V,[],1))./sqrt(nSub);
    b = barwitherr(SEM',M'); 
    b(1). FaceColor = [0 0 0.5];
    b(2).FaceColor = [0.5 0 0]; 
    legend('Congruent','Incongruent'); %ylim([1 4])
    set(gca,'XTickLabels',{'Absent','Present'}); ylabel('Vividness');
    for p = 1:2
    hold on; scatter([p-0.15 p+0.15]+randn(nSub,1)./20,squeeze(V(:,:,p)),40,...
            'filled','MarkerFaceColor','k','MarkerFacealpha',0.5)
    end
   
    % plot vividness per response
    figure; 
    nan_idx = squeeze(any(any(any(isnan(Vpr),2),3),4));
    for p = 1:2
        subplot(2,1,p)
        dat = squeeze(Vpr(:,:,:,p)); dat(nan_idx,:,:) = [];
        SEM = squeeze(std(dat)./sqrt(sum(~nan_idx)));
        M   = squeeze(mean(dat));
        b = barwitherr(SEM,M);
        for i = 1:2
            hold on; scatter(((1:2)+b(i).XOffset)+(randn(sum(~nan_idx),1)./20),...
                squeeze(dat(:,:,i)),30,'k','filled','MarkerFaceAlpha',0.5)
        end
        set(gca,'XTickLabel',{'congruent','incongruent'});
        ylim([1.5 4])
    end

    
end
