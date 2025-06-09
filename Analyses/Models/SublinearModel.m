function [rj_sim,v_sim] = SublinearModel(theta,Vm,Dp,Pres,Cond)

nTrls = length(Pres);
rj_sim = nan(nTrls,1);
v_sim = nan(nTrls,1);

c = 1000; % constant with to add to ensure non negative numbers for power function 

% convert theta to the right parameters
alpha = 1/(1+exp(-theta(1)))*4;
beta  = 1/(1+exp(-theta(2)))*4;
gamma = 1/(1+exp(-theta(3))); % between 0 and 1

for t = 1:nTrls
    I = mvnrnd([Vm 0],[1 0; 0 1]);
    if Cond(t)==1
        if Pres(t)==1
            P = mvnrnd([Dp 0],[1 0; 0 1]);
        elseif Pres(t)==0
            P = mvnrnd([0 0],[1 0; 0 1]);
        end

        RS = (alpha.*P)+(beta.*I);
        RS = (RS+c).^gamma; % sublinear - add constant to prevent negativity            

        rj_sim(t) = RS(1);
        v_sim(t) = RS(1);
    elseif Cond(t)==2
        if Pres(t)==1
            P = mvnrnd([0 Dp],[1 0; 0 1]);
        elseif Pres(t)==0
            P = mvnrnd([0 0],[1 0; 0 1]);
        end

        RS = (alpha.*P)+(beta.*I);
        RS = (RS+c).^gamma; % sublinear - add constant to prevent negativity            

        rj_sim(t) = RS(2);
        v_sim(t) = RS(1);
    end
end