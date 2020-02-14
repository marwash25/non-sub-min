function Av_pinv = pinvupdate(A,A_pinv,v)
%% function to update pseudo inverse of [A | v]  pseudo inverse A_pinv of A
% [A |v] = [A | 0] + v*[0...0 | 1]
% pinv([A | 0]) = [pinv(A) // 0]

[m,n] = size(A);
A0 = [A, zeros(m,1)];
A0_pinv = [A_pinv;zeros(1,m)];
d = [zeros(n,1);1];

%tol = max(size(A))*eps(norm(A));
Av_pinv = OneRankInverseUpdate(A0,A0_pinv,v,d);

end