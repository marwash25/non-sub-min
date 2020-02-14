function [x_primal,x_dual,dual_values,primal_values,gaps,added1,added2,added3] = minimize_submodular_FW_condgradient_nolinesearch(F,F_true,param_F,maxiter,sfm, optimal_primal)
% Compute  min_x .5 * || x ||^2 + f(x) using Frank-Wolfe on the dual with line
% search (i.e., conditional gradient) - NO LINE SEARCH
%
% if sfm=1, perform submodular function minimization (and adapt gaps and
% function values accordingly)
%
%
% INPUT
% F,param_F: submodular function
% maxiter: maximum number of iterations
% sfm: 1 if minimizing submodular function, 0 otherwise
%
% OUTPUT
% x: argmin. If sfm=1, then outputs the subset
% values of cost function
% gaps: certified optimaly gaps at each iteration

if nargin<4, sfm=0; end

if nargout>=8
    tic;
end

p = param_F.p;

FV = F(1:p,param_F);
Fsingletons = zeros(p,1);
FVsingletons = zeros(p,1);

for i=1:p
    Vmi = 1:p; Vmi(i)=[];
    FVsingletons(i) = FV-F(Vmi,param_F);
    Fsingletons(i) = F(i,param_F);
end


% random initialization
direction = rand(p,1);
direction((Fsingletons<=0)) = 1;
direction((FVsingletons>=0)) = 0;
[s,Fvalues,order] = greedy_algo_submodular(direction,F,param_F);
bestvalue = Inf;

for iter =1:maxiter
    
    % order return the descending order of (-s) so the ascending order of
    % s, and Fvalues corresponds to the values of all the possible sublevel 
    % sets except the empty set (i.e., p sets)
    [w,Fvalues,order] = greedy_algo_submodular(-s,F,param_F);
    
    % compute certificates
    if sfm
        % a is the sublevel set with the minimum value
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
        primal_values(iter)  = F_true(Aopt,param_F); %a;
        dual_values(iter) = sum(min(s,0));
        gaps(iter) = a - dual_values(iter) ;
        
        % compute traditional value
        added1(iter) = F(find(s<-1e-15),param_F);
        added2(iter) = F(find(s<+1e-15),param_F);
        
        if primal_values(iter) < optimal_primal
            fprintf('stop because reached required accuracy on SFM %f\n',optimal_primal - primal_values(iter)); 
            break;
        end
    else
        if nargout>=8
            added3(iter) = toc;
        end
        gaps(iter) = -s'*w + s'*s;
        dual_values(iter) =  -.5 * (s'*s);
        ww = pav( -w(flipud(order)) );
        ww(flipud(order)) = ww;
        added1(iter) = w'* ( - s) + .5 * s'*s;
        primal_values(iter) = ww'*w + .5 * ww' * ww;
        added2(iter) = added1(iter) - dual_values(iter);
        gaps(iter) = primal_values(iter) - dual_values(iter);
        x_primal = ww;
    end
    
    
    % no line search
    s = s + (w-s) * 2/(2 + iter);
    
    
end

% output
x_dual = s;
if sfm
    x_primal = Aopt;
end

