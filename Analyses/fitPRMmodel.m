function [params,fVal] = fitPRMmodel(RJ,Vt,Vm,Cond,Pres,Dp,theta,m)

options = optimset;
[params, fVal] = fminsearch(@fitfunc, theta, options);

disp(['dev = ' num2str(fVal)])
disp(['alpha = ' num2str(1/(1+exp(-params(1)))*4')])
disp(['beta = ' num2str(1/(1+exp(-params(2)))*4')])

    function dev = fitfunc(theta)

        % convert theta to alpha and beta
        alpha = 1/(1+exp(-theta(1)))*4;
        beta  = 1/(1+exp(-theta(2)))*4;

        rng(1); % reset rng
        nIter = 100;
        nTrls    = length(RJ);
        rj_sim = nan(nTrls,nIter);
        v_sim  = nan(nTrls,nIter);

        for i = 1:nIter

            if m == 1 % no mixing/null-model

                for t = 1:nTrls
                    if Cond(t) == 1 % congruent
                        if Pres(t)==1
                            P = mvnrnd([Dp 0],[1 0; 0 1]);
                        elseif Pres(t)==0
                            P = mvnrnd([0 0],[1 0; 0 1]);
                        end
                        P = alpha.*P;
                        rj_sim(t,i) = P(1);
                    elseif Cond(t) == 2 % incongruent
                        if Pres(t)==1
                            P = mvnrnd([0 Dp],[1 0; 0 1]);
                        elseif Pres(t)==0
                            P = mvnrnd([0 0],[1 0; 0 1]);
                        end
                        P = alpha.*P;
                        rj_sim(t,i) = P(2);
                    end
                    I = mvnrnd([Vm 0],[1 0; 0 1]);
                    I = beta.*I;
                    v_sim(t,i)  = I(1);
                end


            elseif m == 2 % imagery and perception mixing model

                for t = 1:nTrls
                    I = mvnrnd([Vm 0],[1 0; 0 1]);
                    if Cond(t)==1
                        if Pres(t)==1
                            P = mvnrnd([Dp 0],[1 0; 0 1]);
                        elseif Pres(t)==0
                            P = mvnrnd([0 0],[1 0; 0 1]);
                        end
                        RS = (alpha.*P)+(beta.*I);
                        rj_sim(t,i) = RS(1);
                        v_sim(t,i) = RS(1);
                    elseif Cond(t)==2
                        if Pres(t)==1
                            P = mvnrnd([0 Dp],[1 0; 0 1]);
                        elseif Pres(t)==0
                            P = mvnrnd([0 0],[1 0; 0 1]);
                        end
                        RS = (alpha.*P)+(beta.*I);
                        rj_sim(t,i) = RS(2);
                        v_sim(t,i) = RS(1);
                    end
                end

            end
        end

        % average over iterations to get rid of randomness
        rj_sim = mean(rj_sim,2);
        v_sim  = mean(v_sim,2);

        % calculate the fit
        [b_rj, dev_rj, stats] = glmfit(rj_sim, RJ, 'binomial', 'link', 'logit');
        [b_v, dev_v, stats]   = mnrfit(v_sim, Vt, 'model', 'ordinal', 'link','logit');

        % compute total deviance (-2*LL):
        dev = dev_rj + dev_v;
    end

end