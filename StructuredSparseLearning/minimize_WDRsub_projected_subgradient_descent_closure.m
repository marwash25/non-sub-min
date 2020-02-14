%% Slightly modified version of function from Submodular package by Francis Bach
% https://www.di.ens.fr/~fbach/submodular/, retrieved September 19, 2016
%%
function [x_primal,x_dual,dual_values,primal_values,gaps,added1,added2,time,errValues,alphasF, betasG]  = minimize_WDRsub_projected_subgradient_descent_closure(F,G_next,param_H,maxiter,gap_tol,ErrFcns,compute_alphaF,compute_betaG)
% Compute  min_{x in [0,1]^p} h(x) using projected subgradient descent
% H(S) = F(S) - G(S)
%
% INPUT
% F,param_F: submodular function
% maxiter: maximum number of iterations
%
% OUTPUT
% x: argmin F(S)
% gaps: certified optimaly gaps at each iteration


p = param_H.p;

if nargout>=8
    tic
end

L = param_H.L; % set to 3*max_S |F(S)| if known, and Inf otherwise 

%FV = F(1:p,param_F);
Hsingletons = zeros(p,1);
%FVsingletons = zeros(p,1);

for i=1:p
   %Vmi = 1:p; Vmi(i)=[];
   %FVsingletons(i) = FV-F(Vmi,param_F);
   Hsingletons(i) = F(i) - G_next([],i,[]);
end
% 
if L == Inf
    L = sqrt( sum( Hsingletons.^2 ) ); 
end


w = rand(p,1);

% known in advance! 
w(Hsingletons<=0)=1;
%w(FVsingletons>=0)=0;

R = 2*sqrt(p);

bestvalue = Inf;
s_ave = 0;


for iter =1:maxiter

    [s_F, s_G,Hvalues,order] = greedy_algo_WDRSub(w,F,G_next);
    s = s_F - s_G;
    s_ave = (s_ave * (iter-1) + s)/(iter);
    
    % compute function values
    [a,b] = min(Hvalues);
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
        time(iter) = toc;
    end
    added1(iter)  = s'*w;
    primal_values(iter) = min(a,0);
    dual_values(iter) = sum(min(s_ave,0));
    gaps(iter) = primal_values(iter) - dual_values(iter);
    added2(iter) = added1(iter) - dual_values(iter);
    errValues(iter,:) = ErrFcns(Aopt);
    alphasF(iter) = compute_alphaF(s_F);
    betasG(iter) = compute_betaG(s_G);
    
    if gaps(iter) < gap_tol
        fprintf('stop because reached duality gap on SFM %f\n',gaps(iter)); 
        break; 
    end
    % projected gradient descent
    eta = R/(L*sqrt(iter));
    w = w - eta * s;
    w = min(max(w,0),1);
    
end
x_primal = Aopt;
x_dual = s_ave;


