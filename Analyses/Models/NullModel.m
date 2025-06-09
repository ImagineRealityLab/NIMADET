function [rj_sim,v_sim] = NullModel(theta,Vm,Dp,Pres,Cond)

nTrls = length(Pres);
rj_sim = nan(nTrls,1);
v_sim = nan(nTrls,1);


% convert theta to the right parameters
alpha = 1/(1+exp(-theta(1)))*4; % between 0 and 4
beta  = 1/(1+exp(-theta(2)))*4; % between 0 and 4

for t = 1:nTrls
    if Cond(t)==1
        if Pres(t)==1
            P = mvnrnd([Dp 0],[1 0; 0 1]);
        elseif Pres(t)==0
            P = mvnrnd([0 0],[1 0; 0 1]);
        end
        P = alpha.*P;
        rj_sim(t) = P(1);
    elseif Cond(t)==2
        if Pres(t)==1
            P = mvnrnd([0 Dp],[1 0; 0 1]);
        elseif Pres(t)==0
            P = mvnrnd([0 0],[1 0; 0 1]);
        end
        P = alpha.*P;
        rj_sim(t) = P(2);
    end
    I = mvnrnd([Vm 0],[1 0; 0 1]);
    I = beta.*I;
    v_sim(t)  = I(1);
end
