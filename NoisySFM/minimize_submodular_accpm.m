function [primal_variable,dual_variable, dual_values, primal_values, gaps,added1,added2,added3 ] = minimize_submodular_accpm(F,F_true,param_F,maxiter,gap,optimal_primal)
% minimize submodular function by ACCPM: min u s.t f_L(w)<= u
% this is implementing the cutting plane method from Bach-learning survey
% see section 7.5 and 10.6

if nargin<4
    gap = 1e-6;
end

if nargout>=8
    tic
end

n = param_F.p;
%parallel = param_F.parallel;

% use random initialization of permutation
% w = rand(n,1);
% I'm choosing w = 0, to make sure it satisfies f_L(w) <= 0
w = zeros(n,1);
S = greedy_algo_submodular(w,F,param_F);
u = 0;
bestvalue = Inf;
eta = 1;

for iter =1:maxiter
    
    alpha = 1;
    beta = 1;
    epsilon_newton = 1e-12;
    maxiter_newton = 100;
    [eta, lambda2, iter_newton] = accpm_submodular_step(eta,u,S,alpha,beta,epsilon_newton,maxiter_newton);
    s = beta * ( S*eta);
    w = .5  - .5 * s ./ ( 1 + sqrt(1+s.*s));
    
    [news,Fvalues,order] = greedy_algo_submodular(w,F,param_F);
    
    [a,b] = min(Fvalues);
    u = min(u, a);
    s = S*eta/sum(eta);
    
     % compute Aopt
    if min(a,0) < bestvalue
        if a < 0
            % allow empty set to be optimal
            Aopt = order(1:b);
            bestvalue = a;
        else
            bestvalue = 0;
            Aopt = [];
        end
    end
    
    if nargout>=8
        added3(iter) = toc;
    end
    dual_values(iter) = sum(min(s,0));
    primal_values(iter) = F_true(Aopt,param_F); %u;
    gaps(iter) = u - dual_values(iter);
    S = [ S, news];
    eta = [ eta; 1];
    
    if gaps(iter) < gap
        fprintf('stop because reached duality gap on SFM %f\n',gaps(iter));
        break; 
    end
    if primal_values(iter) < optimal_primal
        fprintf('stop because reached required accuracy on SFM %f\n',optimal_primal - primal_values(iter)); 
        break;
    end
end
primal_variable = Aopt;
dual_variable = s;
added1 = [];
added2 = []; 