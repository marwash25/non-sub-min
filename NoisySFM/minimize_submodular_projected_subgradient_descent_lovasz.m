function [x_primal,x_dual,dual_values,primal_values,gaps,added1,added2,added3]  = minimize_submodular_projected_subgradient_descent_lovasz(F,F_true,param_F,maxiter,optimal_primal)
% Compute  min_{x in [0,1]^p} f(x) using projected subgradient descent
%
% if sfm=1, perform submodular function minimization (and adapt gaps and
% function values accordingly)
%
%
% INPUT
% F,param_F: submodular function
% maxiter: maximum number of iterations
%
% OUTPUT
% x: argmin. If sfm=1, then outputs the subset
% values of cost function
% gaps: certified optimaly gaps at each iteration

if nargout>=8
    tic
end

p = param_F.p;
L = param_F.L; % set to 3*max_S |F(S)| if known, and Inf otherwise 

FV = F(1:p,param_F);
Fsingletons = zeros(p,1);
FVsingletons = zeros(p,1);

for i=1:p
    Vmi = 1:p; Vmi(i)=[];
    FVsingletons(i) = FV-F(Vmi,param_F);
    Fsingletons(i) = F(i,param_F);
end

if L == Inf
    L = sqrt( sum( Fsingletons.^2 ) ); 
end
%B = sqrt( sum( (FVsingletons - Fsingletons).^2 ) ); 

R = 2*sqrt(p);

w = rand(p,1);

% known in advance! 
w(Fsingletons<=0)=1;
w(FVsingletons>=0)=0;

bestvalue = Inf;
s_ave = 0;


for iter =1:maxiter

    [s,Fvalues,order] = greedy_algo_submodular(w,F,param_F);
    s_ave = (s_ave * (iter-1) + s)/(iter);
    
    % compute function values
    [a,b] = min(Fvalues);
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
    added1(iter)  = s'*w;
    primal_values(iter) = F_true(Aopt,param_F); %a;
    dual_values(iter) = sum(min(s_ave,0));
    gaps(iter) = a - dual_values(iter);
    added2(iter) = added1(iter) - dual_values(iter);
    
    if primal_values(iter) < optimal_primal
        fprintf('stop because reached required accuracy on SFM %f\n',optimal_primal - primal_values(iter)); 
        break;
    end
    
    % projected gradient descent
    eta = R/(L*sqrt(iter));
    w = w - eta * s;
    w = min(max(w,0),1);
    
end
x_primal = Aopt;
x_dual = s_ave;


