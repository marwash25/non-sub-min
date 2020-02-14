function [x_CR] = minimize_convRel_range(y,A,B,b,w,lbd)

d = size(A,2);
r = size(B,1);

cvx_begin quiet
    cvx_precision best
    variables x_CR(d) z(r);

    minimize (0.5*square_pos(norm(y - A*x_CR, 2))+lbd*b'*z)
    subject to
        B'*z >= abs(w.*x_CR)
        z>= 0
        sum(z) <= 1
cvx_end
                    
end