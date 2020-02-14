function     [eta, lambda2, iter] = accpm_submodular_step(eta,u,S,alpha,beta,epsilon_newton,maxiter_newton)
% perform Newton steps to find weighted analytic center 
% this is solving the dual maximization problem (see Section 10.6 in Bach-Leaning)

% parameters of Newton step
alpha_backtrack = .05;
beta_backtrack = .5;
n = length(eta);
p = size(S,1);
val = max(mean(S.^2,1));

% tau = 1e-3;
% u = u + tau*p*val; %noisy u

for iter=1:maxiter_newton
    % compute gradient, hessian and functin values
    
    [tempf,tempg,temph]=phi_accpm(beta * (S*eta) );
    fx = beta * sum(log(eta)) + alpha * log(sum(eta)) - u*beta*sum(eta) + sum(tempf); 
    gradient = beta ./ eta + alpha/(sum(eta)) - u*beta + beta * S' * tempg;
    hessianx = -beta * diag(1./eta./eta) - alpha/(sum(eta))^2 + beta*beta * S'*diag(temph)*S;
    
    
    
    ascent =  - ( hessianx - val*1e-14*eye(n) )\ gradient;
    lambda2(iter) =  gradient' * ascent;
    
    
    if lambda2(iter) < 2 * epsilon_newton 
        % converged!
        break;
        
    else
        
        % simlpe backtrackinh line search
        a = 0;
        b = Inf;
        kmax_ls = 30;
        kls = 0;
        alpha_bt = 1;
        while kls < kmax_ls
            neweta = eta + alpha_bt * ascent;
            if any(neweta<1e-12), alpha_bt = alpha_bt * beta_backtrack ;
            else
                falpha = beta * sum(log(neweta)) + alpha * log(sum(neweta)) - u*beta*sum(neweta) + sum(phi_accpm(beta * (S*neweta) ));
                if (falpha < fx + alpha_backtrack * lambda2(iter) * alpha_bt)
                    % alpha is too large, reduce it
                    
                    alpha_bt = alpha_bt * beta_backtrack ;
                    %   if display, fprintf('-'); end
                else
                    
                    break;
                end
            end
            kls = kls + 1;
        end
        eta = eta + alpha_bt * ascent;
    end
end
if lambda2(iter) < 0
    fprintf('negative lambda2!');
    %keyboard; 
end


function [f,g,h]=phi_accpm(s);
temp = sqrt(s.*s+1);
f = s/2 + .5* ( 1 - temp) + .5 * log( (1 +temp)/2);
g = .5 - .5 .* s./ (1+temp);
h = -.5 ./ temp ./  (1+temp);

