function x = LS(S,A,y)
x = zeros(size(A,2),1);
%x(S) = A(:,S)\y;
x(S) = pinv(A(:,S))*y;
end