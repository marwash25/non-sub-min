%% Slightly modified version of function from Submodular package by Francis Bach
% https://www.di.ens.fr/~fbach/submodular/, retrieved September 19, 2016
%%
function [x_F,x_G,Hvalues,I] = greedy_algo_WDRSub(w,F,G_next)
% Given a submodular function F(set,param_F) 
% returns x=gradient of the lovasz extension
% x \in argmax w^Ts , s \in B(F)
% Fvalues=values of F(j1 .. jk), where w_j1>= ...>= w_jk
% I= j1 ... jp (sorted coefficients of w)

n = length(w);
x_F = zeros(n,1); 
x_G = zeros(n,1); 
%x_H = zeros(n,1); 
[~,I] = sort(w,'descend');

Fvalues = zeros(n,1);
Gvalues = zeros(n,1);
%Hvalues = zeros(n,1);

% Fvalues(1) = F(I(1),param_F);
% x(I(1)) = Fvalues(1);
% for i=2:n
%     Fvalues(i) = F(I(1:i),param_F);
%     x(I(i)) = Fvalues(i) - Fvalues(i-1);
% end

Fvalues(1) = F(I(1));
[Gvalues(1),AS_pinv] = G_next([],I(1),[]);
%[Hvalues(1),H_AS_pinv] = H_next([],I(1),[]);

x_F(I(1)) = Fvalues(1);
x_G(I(1)) = Gvalues(1);
%x_H(I(1)) = Hvalues(1);

for i=2:n
    Fvalues(i) = F(I(1:i));
    [Gvalues(i),AS_pinv] = G_next(I(1:i-1),I(i),AS_pinv);
    %[Hvalues(i),H_AS_pinv] = H_next(I(1:i-1),I(i),H_AS_pinv);
    
    x_F(I(i)) = Fvalues(i) - Fvalues(i-1);
    x_G(I(i)) = Gvalues(i) - Gvalues(i-1);
    %x_H(I(i)) = Hvalues(i) - Hvalues(i-1);
end
Hvalues = Fvalues - Gvalues;