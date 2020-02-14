function [HSe,ASe_pinv] = H_update(S,e,AS_pinv,l,A,y,F)
    [x_LS_Se, ASe_pinv] = LSupdate(S,e,A,y,AS_pinv);
    [~,n] = size(A);
    HSe = F([S;e]) - (l(zeros(n,1)) - l(x_LS_Se));
end