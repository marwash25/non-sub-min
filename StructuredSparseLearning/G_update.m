function [GSe,ASe_pinv] = G_update(S,e,AS_pinv,l,A,y)
    [x_LS_Se, ASe_pinv] = LSupdate(S,e,A,y,AS_pinv);
    [~,n] = size(A);
    GSe = (l(zeros(n,1)) - l(x_LS_Se));
end