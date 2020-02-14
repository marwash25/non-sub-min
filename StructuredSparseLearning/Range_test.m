function Range_test(run,d,k,ind_n,sig_noise,corr_level,maxiter)

%pc = parcluster('local');
%parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
%% Parameters

% set seed
rng(run)
%1- Brute force with range 2- Brute force with modified range 2- SFMin with range, 3- SFMin with modified range, 4- Convex Relaxation for range
algos = [1,2,3,4,5]; 
compute_opt = 1;

n_vec = ceil(linspace(0.25,2,10)*d); 
n = n_vec(ind_n);

lbd_vec = logspace(-4,1,15); 
lbd_vec = sort(lbd_vec,'descend');

constant = 1;
gap_tol = 1e-7;

% Generate all possible intervals I
B = zeros(d*(d+1)/2+1,d);
b = zeros(d*(d+1)/2+1,1); % b = F(I) for all I.
r = 1; 
b(1) = 0; % for the empty set 
for i = 1:d
    for j = i:d
       r = r+1;
       B(r,i:j) = 1;
       b(r) = j-i + 1;
    end
end

%% Generate Gaussian/Constant signal with interval support 

x_true = zeros(d,1);

J_start = randsample(d-k+1,1);
J = J_start:J_start+k-1;
if constant
    %xJ = sign(randn(k,1)); 
    xJ = ones(k,1);
else
    xJ = randn(k,1);
    xJ = xJ/norm(xJ,inf); % normalize to have ||x||inf =1;
end
x_true(J) = xJ; 

EstErr = @(xhat) norm(xhat - x_true,2)/norm(x_true);
SuppErr = @(Jhat) numel(setdiff(J,Jhat)) + numel(setdiff(Jhat,J)); 

%for ind_n = 1:numel(n_vec)
    n  = n_vec(ind_n);

    if corr_level>0
        c = corr_level; % ~correlation b/w cols
        L = triu(c + .01*randn(d),1);
        Q = L + L' + eye(d);
        min_eig  = min(eig(Q));
        if min_eig <0 %make sure Q is psd
            Q = (Q - (min_eig-1e-7)*eye(d))/(1 - (min_eig-1e-7));
        end
        R = chol(Q);
        A = 1/sqrt(n) * randn(n, d)*R;
    else
        A = 1/sqrt(n) * randn(n, d);
    end

    Ax_true = A*x_true;

    noise = sig_noise*randn(n,1);
    y = Ax_true + noise;

    l = @(x) 0.5*norm(y - A*x)^2;
    gradl = @(x) A'*(A*x - y);

    param_H.p = d;
    param_H.L = Inf;
    
    G_next = @(S,e,AS_pinv) G_update(S,e,AS_pinv,l,A,y);
     
    ErrFcns = @(S) [SuppErr(S), EstErr(LS(S,A,y))];
    

%%
    if compute_opt
        time = tic;
        % compute G(S) for all possible intervals
        G_intervals = zeros(d*(d+1)/2+1,1);
        F_intervals = zeros(d*(d+1)/2+1,1);
        F_mod_intervals = zeros(d*(d+1)/2+1,1);

        r = 1; %for the empty set
        G_intervals(1) = 0;
        F_intervals(1) = 0;
        F_mod_intervals(1) = 0;

        for i = 1:d 
            Aij_pinv = [];
            I = [];
            for j = i:d
                r = r+1;
                [G_intervals(r),Aij_pinv] = G_update(I,j,Aij_pinv,l,A,y);
                I = [I;j];
                F_intervals(r) = range_mod(I,0);
                F_mod_intervals(r) = range_mod(I,d);
            end
        end
        time_opt = toc(time);
        
        for ind_algo = 1:2
            for ind_lbd = 1:length(lbd_vec)
                lbd = lbd_vec(ind_lbd);
                time2 = tic;
                if ind_algo ==1
                    %compute true minimizer of l(x) + lbd*range(supp(x)) by brute force
                    [H_opt, ind_opt] = min(lbd*F_intervals - G_intervals);    
                    F_opt_values(ind_lbd) = F_intervals(ind_opt);
                    G_opt_values(ind_lbd) = G_intervals(ind_opt);
                else
                    %compute true minimizer of l(x) + lbd*range_mod(supp(x)) by brute force
                    [H_opt, ind_opt] = min(lbd*F_mod_intervals - G_intervals);
                    F_mod_opt_values(ind_lbd) = F_mod_intervals(ind_opt);
                    G_mod_opt_values(ind_lbd) = G_intervals(ind_opt);
                end
                S_opt = find(B(ind_opt,:)~=0);
                time2 = toc(time2);
                time_values(ind_lbd) = time_opt + time2;

                x = LS(S_opt,A,y);
                S_values{ind_lbd} = S_opt;
                x_values{ind_lbd} = x;
                H_values(ind_lbd) = H_opt;
                dual_values(ind_lbd) = sum(min(x ,0));
                SuppErr_values(ind_lbd) = SuppErr(S_opt);
                EstErr_values(ind_lbd) = EstErr(x);
            end
            
        results{ind_algo}.ObjErr = H_values;
        results{ind_algo}.dual = dual_values;
        results{ind_algo}.SuppErr = SuppErr_values;
        results{ind_algo}.EstErr = EstErr_values;
        results{ind_algo}.time = time_values;
        results{ind_algo}.S = S_values;
        results{ind_algo}.x = x_values;
        results{ind_algo}.alphaF = 1;
        results{ind_algo}.betaG = 1;
        results{ind_algo}.nits = r;
        end
        
    end

            
    for ind_algo = 3:length(algos) %run algos sequentially 
        algo = algos(ind_algo);
        fprintf('Running algorithm %d \n',algo)

        parfor ind_lbd = 1:length(lbd_vec) %inefficient. TO DO: fix it!
            lbd = lbd_vec(ind_lbd);

            switch algo
                    
                case 3 
                    % the sol we find should satisfy H_opt <= H(S) <= min{F_opt/alpha_F,F(V)} - beta_G G_opt 
                    % alpha_F = 1/(d-1) here, beta_G = min_{S<=T,i not in T} mu_|S| grad(x*(T))^2_i / nu_|T| grad(x*(S))^2_i 
                    
                    F = @(S) lbd*range_mod(S,0); % range function 
                    %H_next = @(S,e,AS_pinv) H_update(S,e,AS_pinv,l,A,y,F);
                    compute_alphaF = @(s) min(F_opt_values(ind_lbd)/sum(s(results{1}.S{ind_lbd})),1); %alpha_F <= F(OPT)/kappa_t(OPT)
                    compute_betaG = @(s) min(sum(s(results{1}.S{ind_lbd}))/G_opt_values(ind_lbd),1); %beta_G <= kappa_t(OPT)/G(OPT)
                    [S,~,dual,primal,gaps,~,~,time, errors,alphasF, betasG]  = minimize_WDRsub_projected_subgradient_descent_closure(F,G_next,param_H,maxiter,gap_tol, ErrFcns, compute_alphaF,compute_betaG);
                    x = LS(S,A,y);

                    S_values{ind_lbd} = S;
                    x_values{ind_lbd} = x;
                    [H_best,ind_best] = min(primal);
                    H_values(ind_lbd) = H_best;
                    time_values(ind_lbd) = time(ind_best);
                    nits_values(ind_lbd) = ind_best;
                    dual_values(ind_lbd) = min(dual);
                    SuppErr_values(ind_lbd) = min(errors(:,1));
                    EstErr_values(ind_lbd) = min(errors(:,2));
                    alphaF(ind_lbd) = mean(alphasF); %this is the alpha that will show up in the bound
                    betaG(ind_lbd) = mean(betasG); %this is the beta that will show up in the bound
                    
                case 4 
                   %the sol we find should satisfy H_mod_opt <= H_mod(S) <= F_mod_opt - beta_G G_opt 
                   % beta_G = min_{S_opt<=T,i not in T} mu_|S_opt| grad(x*(T))^2_i / nu_|T| grad(x*(S_opt))^2_i

                    F_mod = @(S) lbd*range_mod(S,d); % modified range function 
                    H_mod_next = @(S,e,AS_pinv, param_H) H_update(S,e,AS_pinv,l,A,y,F_mod);
                    compute_alphaF = @(s) min(1,F_mod_opt_values(ind_lbd)/sum(s(results{2}.S{ind_lbd}))); %alpha_F <= F(OPT)/kappa_t(OPT) (should be 1 in this case)
                    compute_betaG = @(s) min(1,sum(s(results{2}.S{ind_lbd}))/G_mod_opt_values(ind_lbd)); %beta_G <= kappa_t(OPT)/G(OPT)
                    [S,~,dual,primal,gaps,~,~,time, errors,alphasF, betasG]  = minimize_WDRsub_projected_subgradient_descent_closure(F_mod,G_next,param_H,maxiter,gap_tol, ErrFcns,compute_alphaF,compute_betaG);
                    x = LS(S,A,y);

                    S_values{ind_lbd} = S;
                    x_values{ind_lbd} = x;
                    [H_best,ind_best] = min(primal);
                    H_values(ind_lbd) = H_best;
                    time_values(ind_lbd) = time(ind_best);
                    nits_values(ind_lbd) = ind_best;
                    dual_values(ind_lbd) = min(dual);
                    SuppErr_values(ind_lbd) = min(errors(:,1));
                    EstErr_values(ind_lbd) = min(errors(:,2));
                    alphaF(ind_lbd) = mean(alphasF); %this is the alpha that will show up in the bound
                    betaG(ind_lbd) = mean(betasG); %this is the beta that will show up in the bound
                    
                case 5
                    F = @(S) lbd*range_mod(S,0); % range function 
                    w = ones(d,1);
                    time = tic;
                    x_CR = minimize_convRel_range(y,A,B,b,w,lbd);
                    
                    %S = find(abs(x_CR) > 1e-5);
                    % round by thresholding
                    [~,~,Hvalues,order] = greedy_algo_WDRSub(abs(x_CR),F,G_next);
                    [H_best,ind_best] = min(Hvalues);
                    if H_best< 0
                        S = order(1:ind_best);
                    else
                        S = [];
                        H_best = 0;
                    end
                    time = toc(time);
                    x = LS(S,A,y);
                    
                    S_values{ind_lbd} = S;
                    x_values{ind_lbd} = x;
                    time_values(ind_lbd) = time;
                    nits_values(ind_lbd) = 1;
                    H_values(ind_lbd) = H_best;
                    dual_values(ind_lbd) = sum(min(x,0));
                    SuppErr_values(ind_lbd) = SuppErr(S);
                    EstErr_values(ind_lbd) = EstErr(x);
                    alphaF(ind_lbd) = 1;
                    betaG(ind_lbd) = 1;
                                      
            end
        end
        
        results{ind_algo}.ObjErr = H_values;
        results{ind_algo}.dual = dual_values;
        results{ind_algo}.SuppErr = SuppErr_values;
        results{ind_algo}.EstErr = EstErr_values;
        results{ind_algo}.time = time_values;
        results{ind_algo}.nits = nits_values;
        results{ind_algo}.S = S_values;
        results{ind_algo}.x = x_values;
        results{ind_algo}.alphaF = alphaF;
        results{ind_algo}.betaG = betaG;


    end
    
%% save results
fprintf('Job completed, saving results :)\n')
corr_str = strrep(num2str(corr_level), '.', '');
sig_str = strrep(num2str(sig_noise), '.', '');
file_name = sprintf('results_range/results%d/RangeTest_corr%s_sigma%s_d%d_k%d_indn%d',run,corr_str,sig_str,d,k,ind_n);
save(file_name,'results')
    

end
 


