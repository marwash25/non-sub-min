function f = submodular_fct_logdetsym(A,param_F)
% logdet of submatrix minus modular function
VA = 1:param_F.p; VA(A) = [];
f = sum(log(diag(chol(param_F.K(A,A)))))*2  + sum(log(diag(chol(param_F.K(VA,VA)))))*2 - param_F.FV;
f = f - sum( param_F.s(A) );
