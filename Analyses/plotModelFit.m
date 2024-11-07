function [AUC, RJ_pred, v_pred] = plotModelFit(RJ,Vt,Vm,Cond,Pres,Dp,theta,m,plotting)

% convert theta to alpha
alpha = 1/(1+exp(-theta))*4;

rng(1); % reset rng
nSS     = 10; % how many sub-samples to equate class sizes?
nTrls    = length(RJ);
nIter   = 100;
RJ_pred = nan(nTrls,nIter);
v_pred  = nan(nTrls,nIter);

%% Get the predicted values
for i = 1:nIter
    if m == 1 % H0: source seperation
        for t = 1:nTrls
            if Cond(t) == 1 % congruent
                if Pres(t)==1
                    P = alpha.*mvnrnd([Dp 0],[1 0; 0 1]);
                elseif Pres(t)==0
                    P = alpha.*mvnrnd([0 0],[1 0; 0 1]);
                end
                RJ_pred(t,i) = P(1)>C;
            elseif Cond(t) == 2 % incongruent
                if Pres(t)==1
                    P = alpha.*mvnrnd([0 Dp],[1 0; 0 1]);
                elseif Pres(t)==0
                    P = alpha.*mvnrnd([0 0],[1 0; 0 1]);
                end
                RJ_pred(t,i) = P(2)>C;
            end
            I = mvnrnd([Vm 0],[0.5 0; 0 0.5]);
            v_pred(t,i)  = I(1);
        end

    elseif m == 2 % H2: complete mixing
        for t = 1:nTrls
            I = mvnrnd([Vm 0],[0.5 0; 0 0.5]);
            if Cond(t)==1
                if Pres(t)==1
                    P = (alpha.*mvnrnd([Dp 0],[1 0; 0 1]))+I;
                elseif Pres(t)==0
                    P = (alpha.*mvnrnd([0 0],[1 0; 0 1]))+I;
                end
                RJ_pred(t,i) = P(1)>C;
                v_pred(t,i) = P(1);
            elseif Cond(t)==2
                if Pres(t)==1
                    P = (alpha.*mvnrnd([0 Dp],[1 0; 0 1]))+I;
                elseif Pres(t)==0
                    P = (alpha.*mvnrnd([0 0],[1 0; 0 1]))+I;
                end
                RJ_pred(t,i) = P(2)>C;
                v_pred(t,i) = P(1);
            end
        end
    end
end

% average over iterations
RJ_pred = mean(RJ_pred,2);
v_pred  = mean(v_pred,2);

% calculate AUC
[x1,y1,~,AUC(1)] = perfcurve(RJ,RJ_pred,1);
[x2,y2,~,AUC(2)] = perfcurve(Vt,v_pred,max(Vt));


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