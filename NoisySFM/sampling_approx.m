function FS = sampling_approx(S, param_F,F_noisy,m, parallel)
% sample F_noisy(S) m times and return the average
    FS_samples = zeros(m,1);
    
    if parallel && (m > 1)
        parfor i = 1:m
            FS_samples(i) = F_noisy(S,param_F);
        end
    else
        for i = 1:m
            FS_samples(i) = F_noisy(S,param_F);
        end
    end
    
    FS = sum(FS_samples)/m;
end