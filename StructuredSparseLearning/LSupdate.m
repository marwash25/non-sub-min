function [x, ASe_pinv] = LSupdate(S,e,A,y,AS_pinv)
x = zeros(size(A,2),1);
ASe_pinv = pinvupdate(A(:,S),AS_pinv,A(:,e));
x([S;e]) = ASe_pinv*y;

end