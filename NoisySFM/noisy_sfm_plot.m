%% plot noisy sfm results
clc
clear all
close all

dataset = 1;
legends = ["MNP","CG-LS","CG-2/(t+2)","PGM-$1/\sqrt{t}$","PGM-polyak","ellipsoid","simplex","ACCPM", "ACCPM-kelley"];
colors = ["k", "r","r--","b","b--","y","c","g","g--"];
algos = [1:5,8,9]; 

switch dataset
    case 1
        name='genrmflong';
        maxiter = 1000;
        d = 575;
    case 2
        name='genrmfwide';
        maxiter = 1500;
        d = 430;
    case 3
        name='twomoons';
        maxiter = 650;
        d = 400;
    case 4
        name='bilmes';
        maxiter = 400;
        d = 800;
    case 6
        name='image';
        maxiter = 400;
        d = 2500;
end
%%
sigma = 0.1;
sig_str = strrep(num2str(sigma), '.', '');

if dataset==3
m_vec_plt =  [1,50,100,150];
else
m_vec_plt =  [1,50,100,300,500,1000];
end

plt = true;
save_fig = false;
run = 1;

%for run = 1:nruns %aveage results over Monte Carlo runs
    for m_ind = 1:length(m_vec_plt)
        m = m_vec_plt(m_ind);
        %fload = sprintf('results_noisySFM/results_laptop/results_dataset%d_samples%d_sig%s',dataset,m,sig_str);
        fload = sprintf('results_noisySFM/results%d/results_dataset%d_samples%d_sig%s.mat',run,dataset,m,sig_str);
        if  exist(fload, 'file')==2
            load(fload)
        end

        for ind_i = 1:length(algos)

            fload = sprintf('results_noisySFM/results%d/results_dataset%d_samples%d_sig%s_algo%d.mat',run,dataset,m,sig_str,ind_i);
            if exist(fload, 'file')==2
                load(fload)
            end

            i = algos(ind_i);   

            total_results{ind_i}.primal = results{ind_i}.primal;
            total_results{ind_i}.dual = results{ind_i}.dual;
            total_results{ind_i}.time = results{ind_i}.time;
           [total_results{ind_i}.primal_best(m_ind), ind_best]= min(results{ind_i}.primal);
            total_results{ind_i}.iter_best(m_ind) = ind_best;
            total_results{ind_i}.time_best(m_ind) = results{ind_i}.time(ind_best);
        end

        if dataset ==1
            m_plt = (m == 50 || m== 1000); 
        elseif dataset ==3
            m_plt = (m == 50 || m== 150);
        end
        
        if plt && m_plt
            figure
            for ind_i = 1:length(algos)
                i = algos(ind_i);   
                % primal values
                semilogy(cummin((1e-12+total_results{ind_i}.primal-optimal_primal)),colors(i),'linewidth',2); hold on
            end

            %l = legend(legends(algos),'Location','NorthEastOutside');
            %set(l,'Interpreter','latex')
            set(gca,'fontsize',25)
            axis([0,1000,0,500])
            if dataset==1
                axis([0,1000,0,500])
            elseif dataset ==3 
                axis([0,600,0,500])
            end
            
            xlabel('iterations');
            ylabel('$H(\hat{S})- H(S^*)$','Interpreter','latex');
            %title(['Objective vs iterations, for m = ',num2str(m)])

            if save_fig           
                fig_name = sprintf("results_noisySFM/primal_iter_dataset%d_samples%d_sig%s",dataset,m,sig_str);
                print(gcf,'-dpdf','-r150',fig_name);
            end
        %%
        if false
            figure
            hold all
            for ind_i = 1:length(algos)
                i = algos(ind_i);  
                % primal values vs time
                plot(total_results{ind_i}.time, cummin(log10(1e-12+total_results{ind_i}.primal-optimal_primal)),colors(i),'linewidth',2);
            end

            l = legend(legends(algos),'Location','NorthEastOutside');
            set(l,'Interpreter','latex')
            set(gca,'fontsize',25)
            xlabel('time');
            %ylabel('F(A)-min(F)');
            ylabel('$\log_{10}(H(\hat{S})- H(S^*))$');
            title(['Objective vs time, for m = ',num2str(m)])

            if save_fig
                fig_name = sprintf("results_noisySFM/primal_time_dataset%d_samples%d_sig%s",dataset,m,sig_str);
                print(gcf,'-dpdf','-r150',fig_name);
            end
         %%   
            figure
                hold all
            for ind_i = 1:length(algos)
                i = algos(ind_i);  
                % dual values 
                plot(cummin(log10(1e-12-total_results{ind_i}.dual+optimal_primal)),colors(i),'linewidth',2); hold on;
            end

            l = legend(legends(algos),'Location','NorthEastOutside');
            set(l,'Interpreter','latex')
            set(gca,'fontsize',25)
            xlabel('iterations');
            ylabel('log_{10}(min(F)-s_-(V))');
            title(['dual vs iterations, for m = ',num2str(m)])

            if save_fig
                fig_name = sprintf("results_noisySFM/dual_iter_dataset%d_samples%d_sig%s",dataset,m,sig_str);
                print(gcf,'-dpdf','-r150',fig_name);
            end
        %%    
            figure
                hold all
            for ind_i = 1:length(algos)
                i = algos(ind_i);  
                % dual values vs time
                plot(total_results{ind_i}.time, cummin(log10(1e-12-total_results{ind_i}.dual+optimal_primal)),colors(i),'linewidth',2); hold on;
            end

            l=legend(legends(algos),'Location','NorthEastOutside');
            set(l,'Interpreter','latex')
            set(gca,'fontsize',25)
            xlabel('time');
            ylabel('log_{10}(min(F)-s_-(V))');
            title(['dual vs time, for m = ',num2str(m)])

            if save_fig
                fig_name = sprintf("results_noisySFM/dual_time_dataset%d_samples%d_sig%s",dataset,m,sig_str);
                print(gcf,'-dpdf','-r150',fig_name);
            end
            
        end
        end
    end
%end
%%
figure
for ind_i = 1:length(algos)
    i = algos(ind_i);   
    % best function value achieved vs number of samples
    semilogy(m_vec_plt, cummin((1e-12 + total_results{ind_i}.primal_best-optimal_primal)),colors(i),'marker','*','linewidth',2); hold on
end

if dataset ==1
    %l = legend(legends(algos),'Location','NorthEastOutside');
    l = legend(legends(algos));
    set(l,'Interpreter','latex')
end

set(gca,'fontsize',25)

if dataset==1
    axis([0,m_vec_plt(end),0,500])
elseif dataset ==3 
    axis([0,m_vec_plt(end),0,500])
end

xlabel('number of samples m');
%ylabel('log_{10}(F(A)-min(F))');
ylabel('$H(\hat{S})- H(S^*)$','Interpreter', 'latex');

if save_fig
    fig_name = sprintf("results_noisySFM/primal_samples_dataset%d_sig%s",dataset,sig_str);
    print(gcf,'-dpdf','-r150',fig_name);
end

if false

%%
figure
hold all
for ind_i = 1:length(algos)
    i = algos(ind_i);   
    % time to reach best solution vs number of samples
    plot(m_vec_plt, total_results{ind_i}.time_best,colors(i),'marker','*','linewidth',2); 
end

l = legend(legends(algos),'Location','NorthEastOutside');
set(l,'Interpreter','latex')
set(gca,'fontsize',25)
xlabel('number of samples');
ylabel('time (sec)');

if save_fig
    fig_name = sprintf("results_noisySFM/time_samples_dataset%d_sig%s",dataset,sig_str);
    print(gcf,'-dpdf','-r150',fig_name);
end

%%
figure
hold all
for ind_i = 1:length(algos)
    i = algos(ind_i);   
    % number of iterations to reach best solution vs number of samples
    plot(m_vec_plt, total_results{ind_i}.iter_best,colors(i),'marker','*','linewidth',2); 
end

l = legend(legends(algos),'Location','NorthEastOutside');
set(l,'Interpreter','latex')
set(gca,'fontsize',25)
xlabel('number of samples');
ylabel('iterations');

if save_fig
    fig_name = sprintf("results_noisySFM/iter_samples_dataset%d_sig%s",dataset,sig_str);
    print(gcf,'-dpdf','-r150',fig_name);
end
end