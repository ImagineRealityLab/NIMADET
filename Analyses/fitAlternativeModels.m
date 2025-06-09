function [params,fVal,output] = fitAlternativeModels(RJ,Vt,Vm,Cond,Pres,Dp,theta,m)
% 1 = Null model
% 2 = Additive model
% 3 = Multiplicative model
% 4 = Sublinear model 
% 5 = Supralinear model

options = optimset;
%options.Display = 'iter';
[params, fVal,~,output] = fminsearch(@fitfunc, theta, options);

    function dev = fitfunc(theta)       

        rng(1); % reset rng
        nIter = 100;
        nTrls    = length(RJ);
        rj_sim = nan(nTrls,nIter);
        v_sim  = nan(nTrls,nIter);

        for i = 1:nIter

            if m == 1 % Null model 
                [rj_sim(:,i),v_sim(:,i)] = NullModel(theta,Vm,Dp,Pres,Cond);

            elseif m == 2 % Additive model                 
                [rj_sim(:,i),v_sim(:,i)] = AdditiveModel(theta,Vm,Dp,Pres,Cond);                 

            elseif m == 3 % Multiplicative
                [rj_sim(:,i),v_sim(:,i)] = MultiplicativeModel(theta,Vm,Dp,Pres,Cond);                    

            elseif m == 4 % Sublinear model
                [rj_sim(:,i),v_sim(:,i)] = SublinearModel(theta,Vm,Dp,Pres,Cond);
            end
        end

        % average over iterations to get rid of randomness
        rj_sim = mean(rj_sim,2);
        v_sim  = mean(v_sim,2);

        % calculate the fit
        warning('off')        
        [b_rj, dev_rj, stats] = glmfit(rj_sim, RJ, 'binomial', 'link', 'logit');
        [b_v, dev_v, stats]   = mnrfit(v_sim, Vt, 'model', 'ordinal', 'link','logit');
        warning('on')

        % compute total deviance (-2*LL):
        dev = dev_rj + dev_v;

        %display(theta);
    end

end