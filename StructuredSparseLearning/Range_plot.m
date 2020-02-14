clc
clear all
close all

legends = ["OPT-Range","OPT-ModRange", "PGM-Range", "PGM-ModRange","CR-Range"];

lineTypes = ["k*-","k--o","r-x","b-+","g-s"];

corr_level = 0;
corr_str = strrep(num2str(corr_level), '.', '');
sig_noise = 0.01;
sig_str = strrep(num2str(sig_noise), '.', '');

d = 250;
k = 20;
algos = [1,2,3,4,5];
n_vec = ceil(linspace(0.25,2,10)*d); 
lbd_vec = logspace(-4,1,15); 
lbd_vec = sort(lbd_vec,'descend');

save_fig = false;
plt = true;

runs = [1:5];

for ind_n = 1:numel(n_vec)
    n  = n_vec(ind_n);
    for ind_run = 1:numel(runs)
        run = runs(ind_run);
    
        fload = sprintf('results_range/results%d/RangeTest_corr%s_sigma%s_d%d_k%d_indn%d.mat',run,corr_str,sig_str,d,k,ind_n);
        if exist(fload, 'file')==2
            load(fload)

             for ind_i = 1:length(algos)

                total_results{ind_i}.ObjErr(ind_run,:,ind_n) = (results{ind_i}.ObjErr);
                %total_results{ind_i}.dual(ind_n) = max(results{ind_i}.dual);
               %[ total_results{ind_i}.SuppErr(ind_run,ind_n), total_results{ind_i}.bestlbd(ind_run,ind_n)] = min(results{ind_i}.SuppErr);
                total_results{ind_i}.SuppErr(ind_run,ind_n) = min(results{ind_i}.SuppErr);
                %total_results{ind_i}.EstErr(ind_run, ind_n) = results{ind_i}.EstErr(total_results{ind_i}.bestlbd(ind_run,ind_n));
               [ total_results{ind_i}.EstErr(ind_run, ind_n), total_results{ind_i}.bestlbd(ind_run,ind_n)]  = min(results{ind_i}.EstErr);
                total_results{ind_i}.alphaF(ind_run,:, ind_n) = results{ind_i}.alphaF;
                total_results{ind_i}.betaG(ind_run, :, ind_n) = results{ind_i}.betaG;
                total_results{ind_i}.time(ind_run, ind_n) = mean(results{ind_i}.time); 
                total_results{ind_i}.nits(ind_run, ind_n) = mean(results{ind_i}.nits);
                total_results{ind_i}.S = results{ind_i}.S;
                %total_results{ind_i}.x = results{ind_i}.x;

             end
             
        else
            fprintf('File does not exist!')
        end
    end
    
    if plt &&  (ind_n == 2 || ind_n == 6)      
        figure
        
        for ind_i = 1:length(algos) %range mod is minimizing a different H
            i = algos(ind_i); 
            semilogx(lbd_vec,mean(total_results{ind_i}.ObjErr(:,:,ind_n),1),lineTypes(i),'markersize',10,'linewidth',2); hold on
        end
%         semilogx(lbd_vec,mean(total_results{1}.bd(:,:,ind_n),1),'--', 'Color',colors(1),'linewidth',2); hold on
%         semilogx(lbd_vec,mean(total_results{2}.bd(:,:,ind_n),1),'--', 'Color',colors(2),'linewidth',2); hold on
        
        xlabel('\lambda (logscale)');
        ylabel('H(S)');
        l = legend(legends(algos),'Location','best');
        %set(l,'Interpreter','latex')
        set(gca,'fontsize',25)
        %title(['Objective vs regularization parameter, for n = ',num2str(n)])
        
        if save_fig &&  (ind_n == 2 || ind_n == 5 || ind_n == 6)      
            fig_name = sprintf("results_range/Objective_lambda_corr%s_sigma%s_d%d_k%d_n%d",corr_str,sig_str,d,k,n);
            print(gcf,'-dpdf','-r150',fig_name);
        end


        figure
         
        %for ind_i = 3:4
            i = algos(ind_i); 
            semilogx(lbd_vec, mean(total_results{3}.alphaF(:,:,ind_n),1),'r--x','markersize',10,'linewidth',2); hold on
            semilogx(lbd_vec, mean(total_results{3}.betaG(:,:,ind_n),1),'r--d','markersize',10,'linewidth',2); hold on
            semilogx(lbd_vec, mean(total_results{4}.alphaF(:,:,ind_n),1),'b--+','markersize',10,'linewidth',2); hold on
            semilogx(lbd_vec, mean(total_results{4}.betaG(:,:,ind_n),1),'b--d','markersize',10,'linewidth',2); hold on
        %end
        ylim([0,1]);
        xlabel('\lambda (logscale)');
        ylabel('parameters');
        
        l = legend(["$\alpha_{T}$-Range", "$\beta_{T}$-Range","$\alpha_{T}$-ModRange", "$\beta_{T}$-ModRange"],'Location','best');
        set(l,'Interpreter','latex')
        set(gca,'fontsize',25)
        %title(['parameters, for n = ',num2str(n)])
        
        if save_fig &&  (ind_n == 2 || ind_n == 5 || ind_n == 6)      
            fig_name = sprintf("results_range/parameters_lambda_corr%s_sigma%s_d%d_k%d_n%d",corr_str,sig_str,d,k,n);
            print(gcf,'-dpdf','-r150',fig_name);
        end
    
    end
   
         
end
%%
figure

for ind_i = 1:length(algos) 
    i = algos(ind_i);   
    plot(n_vec/d, mean(total_results{ind_i}.SuppErr,1),lineTypes(i),'markersize',10,'linewidth',2); hold on
end
xlabel('n/d');
ylabel('Support Error');
xlim([n_vec(1)/d, n_vec(end)/d])
l = legend(legends(algos),'Location','best');
set(l,'Interpreter','latex')
set(gca,'fontsize',25)

if save_fig           
    fig_name = sprintf("results_range/SuppErr_samples_corr%s_sigma%s_d%d_k%d",corr_str,sig_str,d,k);
    print(gcf,'-dpdf','-r150',fig_name);
end


figure

for ind_i = 1:length(algos) 
    i = algos(ind_i);   
    semilogy(n_vec/d, mean(total_results{ind_i}.EstErr,1),lineTypes(i),'markersize',10,'linewidth',2); hold on
end
xlabel('n/d');
ylabel('Estimation Error');


l = legend(legends(algos),'Location','best');
set(l,'Interpreter','latex')
set(gca,'fontsize',25)


if save_fig           
    fig_name = sprintf("results_range/EstErr_samples_corr%s_sigma%s_d%d_k%d",corr_str,sig_str,d,k);
    print(gcf,'-dpdf','-r150',fig_name);
end

