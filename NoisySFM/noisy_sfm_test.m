function noisy_sfm_test(run,dataset,ind_algo,m,sigma,maxiter_in)
%     numCores = feature('numcores')
%     p = parpool(numCores);
    optimal_primal = -Inf;

    switch dataset
        case 1
            name='genrmflong';
            maxiter = 1000;
            d = 575;
            %optimal_primal= -366;
        case 2
            name='genrmfwide';
            maxiter = 1500;
            d = 430;
        case 3
            name='twomoons';
            maxiter = 650;
            d = 400;
            %optimal_primal= -37.8337;
        case 4
            name='bilmes';
            maxiter = 400;
            d = 800;
        case 6
            name='image';
            maxiter = 400;
            d = 2500;
    end
   
    seed = run;
    algos = [1:5,8,9]; 
    gap_tol = 1e-8;
    maxiter = min(maxiter,maxiter_in); %choose min between preset maxiter and input maxiter.
%     m_vec = floor(linspace(1,d/4,10).^2);
%     m = m_vec(ind_m);
    
    fprintf('Running Noisy SFM on dataset %s, with m = %d samples, and sigma = %2.5f noise variance \n\n',name,m,sigma);
    %% compute opt of true function
    if optimal_primal == -Inf

    [true_solution,true_certificate,true_dual,true_primal,true_gap,true_time] = launch_sfm(dataset,0,maxiter,seed,gap_tol,optimal_primal,1,sigma);
    
    optimal_dual = max(true_dual);
    optimal_primal = min(true_primal)
    end
    %% run different algorithms, given access to noisy oracle. Stop when maxiter is reached or primal - optimal_primal <= tol or gap <= gap_tol
    tol = 1e-4;

    %parfor ind_algo = 1:length(algos)
        i = algos(ind_algo);
        fprintf('Running algorithm %d \n',i)
        [solution,certificate,dual,primal,gap,time] = launch_sfm(dataset,i,maxiter,seed,gap_tol,optimal_primal + tol,m,sigma);
        results{ind_algo}.primal = primal;
        results{ind_algo}.dual = dual;
        results{ind_algo}.gap = gap;
        results{ind_algo}.solution = solution;
        results{ind_algo}.certificate  = certificate;
        results{ind_algo}.time = time;
    %end
    
   %% save results
   fprintf('Job completed, saving results :)\n')
   sig_str = strrep(num2str(sigma), '.', '');
   file_name = sprintf('results_noisySFM/results%d/results_dataset%d_samples%d_sig%s_algo%d',run,dataset,m,sig_str,ind_algo);
   save(file_name)

end