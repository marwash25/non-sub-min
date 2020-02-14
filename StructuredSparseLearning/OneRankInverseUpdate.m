% % % Navvab Kashiri (2020). Moore-Penrose Pseudo-Inverse Rank-1 Update 
% % % (https://www.mathworks.com/matlabcentral/fileexchange/61115-moore-penrose-pseudo-inverse-rank-1-update), MATLAB Central File Exchange. Retrieved April 21, 2019.
% % % A function for rank-1 update for the Moore-Penrose pseudo-inverse of real valued matrices
% % % The Matrix Cookbook [ http://matrixcookbook.com ]
% % % 3.2.7 Rank-1 update of Moore-Penrose Inverse
% % % When the matrix A is updated by the product of two vectors c,d: A+cd'
% % % The inputs are :  A, Original Matrix
% % %                   A_pinv, Pseudo Inverse of Matrix A 
% % %                   c, the first input vector
% % %                   d, the second input vector
% % % and the last two optional inputs are: 
% % % A tolerance for checking zero conditions, to conver numerical errors;
% % % A 0/1 flag to print or not print the condition of two input  vectors 
% % % with respect to the matrix A, by setting it to 1 and 0, respectively.
function A_pinv_New = OneRankInverseUpdate(A,A_pinv,c,d,Zero_tol,Case_Print_Flag)
%% Input Checks 
% check the number of inputs
narginchk(4, 6)
% check the dimension of inputs
[n,m] = size(A);
if ~isequal(size(A_pinv),[m,n])
    error('The input pseudo inverse matrix dimension mismatch');
end
if ~isequal(size(c),[n 1])
    error('The first input vector dimension mismatch');
end
if ~isequal(size(d),[m 1])
    error('The second input vector dimension mismatch');
end
% set a small value instead of zero, for avoiding numverical issues
if ~exist('Zero_tol','var')
    Zero_tol = max(size(A))*eps(norm(A));
else 
    if Zero_tol>=1 || Zero_tol<=0
        error('The zero tolerance value should be a very small positive value, and must be between 0 and 1');
    end
end
% set the case number printing status
if ~exist('Case_Print_Flag','var')
    Case_Print_Flag = 0;
else 
    if ~ismember(Case_Print_Flag,[0,1]) 
        error('The last input can only be 0 or 1, i.e. print or hide case number, respectively');
    end
end
%% Initial Calculations
V = A_pinv*c;
b = 1 + d'*V;
N = A_pinv'*d;
W = c - A*V; %(eye(n) - A*A_pinv)*c;
M = d - A'*N; %(eye(m) - (A_pinv*A)')*d;
% squared norm of the two abovesaid vectors
w_snorm = norm(W,2)^2;  
m_snorm = norm(M,2)^2;
%% Computation of the update term 
if w_snorm>=Zero_tol && m_snorm>=Zero_tol
    if Case_Print_Flag == 1
        disp('case 1');
    end
    G = (-1/w_snorm) * V * W' - (1/m_snorm) * M * N' + (b/(m_snorm*w_snorm)) * M * W';
elseif w_snorm<Zero_tol && m_snorm>=Zero_tol && abs(b)<Zero_tol
    if Case_Print_Flag == 1
        disp('case 2'); 
    end    
    v_snorm = norm(V,2)^2;    
    G = (-1/v_snorm) * V * (V'* A_pinv) -(1/m_snorm) * M * N';
elseif w_snorm<Zero_tol && abs(b)>Zero_tol
    if Case_Print_Flag == 1
        disp('case 3');
    end    
    v_snorm = norm(V,2)^2;    
    VtA_pinv = V' * A_pinv;
    G = (1/b) * M * VtA_pinv - (b/(v_snorm * m_snorm + b^2)) *...
        ((v_snorm / b) * M + V) * ((m_snorm / b) * VtA_pinv + N');        
elseif m_snorm<Zero_tol && w_snorm>=Zero_tol && abs(b)<Zero_tol
    if Case_Print_Flag == 1
        disp('case 4');
    end    
    n_snorm = norm(N,2)^2;
    G = (-1/n_snorm) * (A_pinv * N) * N' -(1/ w_snorm) * V * W';
elseif m_snorm<Zero_tol && abs(b)>Zero_tol
    if Case_Print_Flag == 1
        disp('case 5');        
    end    
    n_snorm = norm(N,2)^2;    
    A_pinvN = A_pinv * N;
    G = (1/b) * A_pinvN * W' - (b/(n_snorm * w_snorm + b^2)) *...
        ((w_snorm/b) * A_pinvN + V) * ((n_snorm/b) * W + N)';
elseif m_snorm<Zero_tol && w_snorm<Zero_tol && abs(b)<Zero_tol
    if Case_Print_Flag == 1
        disp('case 6');
    end    
    v_snorm = norm(V,2)^2;
    n_snorm = norm(N,2)^2;
    A_pinvN = A_pinv * N;
    G = (-1/v_snorm) * V * (V' * A_pinv) - (1/n_snorm) * A_pinvN * N' + ...
        (V' * A_pinvN) / (v_snorm * n_snorm) * (V * N');
end
%% Computation of the New Pseudo Inverse
A_pinv_New = A_pinv + G;
end

%%
% Copyright (c) 2017, Navvab Kashiri
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.