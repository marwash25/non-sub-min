 function [x,Fvalues,I] = greedy_algo_submodular(w,F,param_F)
% Given a submodular function F(set,param_F) 
% returns x=gradient of the lovasz extension
% x \in argmax w^Ts , s \in B(F)
% Fvalues=values of F(j1 .. jk), where w_j1>= ...>= w_jk
% I= j1 ... jp (sorted coefficients of w)

n = length(w);
x = zeros(n,1); [temp,I] = sort(w,'descend');
Fvalues = zeros(n,1);

if param_F.parallel
    parfor i = 1:n
       Fvalues(i) = F(I(1:i),param_F);
    end
    x(I(1)) = Fvalues(1);
    %x(I(2:end)) = diff(Fvalues);
    for i=2:n
        x(I(i)) = Fvalues(i) - Fvalues(i-1);
    end
else 
    Fvalues(1) = F(I(1),param_F);
    x(I(1)) = Fvalues(1);
    for i=2:n
        Fvalues(i) = F(I(1:i),param_F);
        x(I(i)) = Fvalues(i) - Fvalues(i-1);
    end
end
