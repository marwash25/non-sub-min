function f = submodular_fct_quad_cut(A,param)
% (semi)-efficient formulation for cut functions
if isempty(A), f=0; else
     A = [ 1 A(:)'+1 ];
    x = sparse(zeros(size(param.W,1),1)); x(A)=1;
    f = x' *  param.W * ( 1 - x ) +  (1-x)' * ( param.W * x ) ;
    x = sparse(zeros(size(param.W,1),1)); x(1)=1;
    f = f - x' *  param.W * ( 1 - x ) -   (1-x)' * ( param.W * x ) ;
    
    
end

