function [BIC, AUC, RJ_pred, v_int, RJ_sim] = modelPredictions(RJ,Vt,Vm,Cond,Pres,Dp,theta,m,plotting)

% convert theta to alpha and beta
alpha = 1/(1+exp(-theta(1)))*4;
beta  = 1/(1+exp(-theta(2)))*4;

rng(1); % reset rng
nTrls    = length(RJ);
nIter   = 100;
RJ_sim = nan(nTrls,nIter);
v_sim  = nan(nTrls,nIter);

%% Get the predicted values
for i = 1:nIter
    if m == 1 % H0: source seperation
        for t = 1:nTrls
            if Cond(t) == 1 % congruent
                if Pres(t)==1
                    P = mvnrnd([Dp 0],[1 0; 0 1]);
                elseif Pres(t)==0
                    P = mvnrnd([0 0],[1 0; 0 1]);
                end
                P = alpha.*P;
                RJ_sim(t,i) = P(1);
            elseif Cond(t) == 2 % incongruent
                if Pres(t)==1
                    P = mvnrnd([0 Dp],[1 0; 0 1]);
                elseif Pres(t)==0
                    P = mvnrnd([0 0],[1 0; 0 1]);
                end
                P = alpha.*P;
                RJ_sim(t,i) = P(2);
            end
            I = mvnrnd([Vm 0],[1 0; 0 1]);
            I = beta.*I;
            v_sim(t,i)  = I(1);
        end

    elseif m == 2 % H2: complete mixing
        for t = 1:nTrls
            I = mvnrnd([Vm 0],[1 0; 0 1]);
            if Cond(t)==1
                if Pres(t)==1
                    P = mvnrnd([Dp 0],[1 0; 0 1]);
                elseif Pres(t)==0
                    P = mvnrnd([0 0],[1 0; 0 1]);
                end
                RS = (alpha.*P)+(beta.*I);
                RJ_sim(t,i) = RS(1);
                v_sim(t,i) = RS(1);
            elseif Cond(t)==2
                if Pres(t)==1
                    P = mvnrnd([0 Dp],[1 0; 0 1]);
                elseif Pres(t)==0
                    P = mvnrnd([0 0],[1 0; 0 1]);
                end
                RS = (alpha.*P)+(beta.*I);
                RJ_sim(t,i) = RS(2);
                v_sim(t,i) = RS(1);
            end
        end
    end
end

% average over iterations
RJ_sim = mean(RJ_sim,2);
v_sim  = mean(v_sim,2);

% calculate the fit
[b_rj, dev_rj, stats] = glmfit(RJ_sim, RJ, 'binomial', 'Link', 'logit');
[b_v, dev_v, stats]   = mnrfit(v_sim, Vt, 'model', 'ordinal', 'link','logit');
dev = dev_rj + dev_v;

% generate predictions by passing through model
RJ_pred = glmval(b_rj,RJ_sim,'logit');
v_pred  = mnrval(b_v,v_sim,'model','ordinal','link','logit');
v_int   = sum(repmat([1:max(Vt)],size(v_pred,1),1).*v_pred,2);

% compute the BIC
if m == 1
    k = 2;
elseif m == 2
    k = 2;
end
BIC = k * log(length(RJ_pred)) - 2 * (-0.5*dev);

% calculate AUC
[x1,y1,~,AUC(1)] = perfcurve(RJ,RJ_pred,1);
[x2,y2,~,AUC(2)] = perfcurve(Vt,v_int,max(Vt));


%% Plot the fit
if plotting
figure;
subplot(1,2,1); % ROC RJ
plot(x1,y1); xlabel('False positive rate'); ylabel('True positive rate');
hold on; plot(0:1,0:1,'k--'); title(sprintf('Reality judgements AUC: %.3f',AUC(1)))

subplot(1,2,2); % ROC viv
plot(x2,y2); xlabel('False positive rate'); ylabel('True positive rate');
hold on; plot(0:1,0:1,'k--'); title(sprintf('Vividness ratings AUC: %.3f',AUC(2)))
end
end