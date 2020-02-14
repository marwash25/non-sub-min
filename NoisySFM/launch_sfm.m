function  [solution,certificate,dual,primal,gap,time] = launch_sfm(dataset,method,maxiter,seed,gap,optimal_primal,m,sigma) 

% launch sfm for some dataset and some method
if nargin<5, gap = 1e-8; end
% load data
[F_true,param_F] = load_data_submodular(dataset,seed);
param_F.L = Inf; %don't know what is max |F(S)|
param_F.parallel = (m > 1);%run function evaluations in parallel inside greedy

% noisy oracle
mu = 1;%mean of the noise
F_noisy = @(S, param_F) (mu + sigma*randn(1))*F_true(S, param_F);

% approximate oracle obtained by sampling
parallel = 0; %making greedy parallel instead
F = @(S, param_F) sampling_approx(S, param_F,F_noisy,m,parallel);


switch method
    
    case 0 % compute the true opt 
        [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_FW_minnormpoint(F_true,F_true,param_F,maxiter,1,gap,optimal_primal);
    
    case 1 % MIN NORM POINT 
        [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_FW_minnormpoint(F,F_true,param_F,maxiter,1,gap,optimal_primal);
        
    case 2 % conditional gradient - line search
        [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_FW_condgradient(F,F_true,param_F,maxiter,1,optimal_primal);
        
    case 3 % conditional gradient - no line search (fixed schedule)
        [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_FW_condgradient_nolinesearch(F,F_true,param_F,maxiter,1,optimal_primal);
        
    case 4 % Subgradient descent - no line search (fixed schedule)
        [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_projected_subgradient_descent_lovasz(F,F_true,param_F,maxiter,optimal_primal);
        
    case 5 % Subgradient descent - Polyak's rule
        [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_projected_subgradient_descent_lovasz_polyak(F,F_true,param_F,maxiter,optimal_primal);
        
%     case 6 % Ellipsoid
%         [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_ellipsoid(F,param_F,maxiter);
%         
%     case 7 % Simplex
%         [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_simplex(F,param_F,maxiter);
%         
    case 8 % ACCPM
        [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_accpm(F,F_true,param_F,maxiter,gap,optimal_primal);

    case 9 % ACCPM - weighted
        [solution,certificate,dual,primal,gap,~,~,time] = minimize_submodular_accpm_weighted(F,F_true,param_F,maxiter,gap,optimal_primal);
       
end


