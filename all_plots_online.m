%% Plot function 
close all;
clc;

n_measurements = 5; %number of observations 
T = 10; %number of time steps 

n_part = 100; %number of particles 250
n_draw = 100; %number of draws for the MC MSE estimator 100

pop_size = 50; %population size of the GA algorithm 
max_gen = 25; %maximum number of generations

n_draw_comp = 40; %number of draws to realize 

try
    load(sprintf('online_computations_matlab_%d_particles_%d_draws',n_part,n_draw));
    n_draw_done = length(find(gain_GA));
    gain_GA = [gain_GA(1:n_draw_done),zeros(1,n_draw_comp)];
    gain_GA_online = [gain_GA_online(1:n_draw_done),zeros(1,n_draw_comp)];
catch
    mse_reg = 0; 
    mse_GA = 0; 
    mse_GA_online = 0; 
    n_draw_done = 0;
    gain_GA = zeros(1,n_draw_comp);
    gain_GA_online = zeros(1,n_draw_comp);
end

mse_reg_comp = 0; 
mse_GA_comp = 0; 
mse_GA_online_comp = 0; 

%% 1. Comparison of the performances of the methods by running n_draw_comp draws
display(['Comparison of the performances of the methods by running draw number ',num2str(n_draw_done+1,'%.0f'),' to draw number ',num2str(n_draw_done + n_draw_comp,'%.0f')]);
t_start_comparison = tic;

meas_reg = round(linspace(0,T,n_measurements));
display('GA computations started');
t_start_GA = tic;

%[meas_GA,~,avgCostHist_GA ,minCostHist_GA] = genetical_algo(n_measurements,T,pop_size,max_gen,2*n_part,2*n_draw);
meas_GA = [1,2,3,6,7];
t_elapsed_GA = toc(t_start_GA);
display(['GA computations completed, time elapsed = ',num2str(t_elapsed_GA,'%.0f sec')]);

measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;

measurements_GA = zeros(1,T+1); 
measurements_GA(meas_GA+1) = 1;


for j = n_draw_done+1:n_draw_done+n_draw_comp
x0 = initialization(1,0,0);
x = model(T,x0,0); %True state vector we want to reconstruct using online OIMPF
y = measurements(x,T,0);

% Computation of the measurement times reg, GA, RT 
display(['Computation of the measurement times GA online, draw number ',num2str(j,'%.0f')]);

display('GA online computations started');
t_start_GA_online = tic;
meas_GA_online = online_optimob(y,n_measurements,T,pop_size,max_gen,n_part,n_draw);
t_elapsed_GA_online = toc(t_start_GA_online);
display(['GA online computations completed, time elapsed = ',num2str(t_elapsed_GA_online,'%.0f sec')]);

n_part_fine = 100000; 

measurements_GA_online = zeros(1,T+1); 
measurements_GA_online(meas_GA_online+1) = 1;

part = initialization(n_part_fine,0,0);
tau_j_reg = particle_filter(y,measurements_reg,T,part,0);
tau_j_GA = particle_filter(y,measurements_GA,T,part,0);
tau_j_GA_online = particle_filter(y,measurements_GA_online,T,part,0);

err_reg = mean((objective(x)-objective(tau_j_reg)).^2,2);
err_GA =  mean((objective(x)-objective(tau_j_GA)).^2,2);
err_GA_online =  mean((objective(x)-objective(tau_j_GA_online)).^2,2);

mse_reg_comp = mse_reg_comp + err_reg;
mse_GA_comp = mse_GA_comp + err_GA; 
mse_GA_online_comp = mse_GA_online_comp + err_GA_online;
        
gain_GA(j) = (err_reg - err_GA)/err_reg;
gain_GA_online(j) = (err_reg - err_GA_online)/err_reg;

save(sprintf('online_computations_matlab_%d_particles_%d_draws',n_part,n_draw),'gain_GA','gain_GA_online')
end

mse_reg = mse_reg*n_draw_done + mse_reg_comp;
mse_GA = mse_GA*n_draw_done + mse_GA_comp;
mse_GA_online = mse_GA_online*n_draw_done + mse_GA_online_comp ;

n_draw_done = n_draw_done + n_draw_comp; %update draw done 

mse_reg = mse_reg/n_draw_done;
mse_GA = mse_reg/n_draw_done;
mse_GA_online = mse_GA_online/n_draw_done;

t_elapsed_comparison = toc(t_start_comparison);

%% 2. Results and save
display(['Computations of ',num2str(n_draw_comp, '%.0f'),'draws, time elapsed = ',num2str(t_elapsed_comparison,'%.0f sec')]);

positive_gain_GA = zeros(n_draw_done,1);
positive_gain_GA(gain_GA>=0) = 1/n_draw_done;
positive_gain_GA = sum(positive_gain_GA);

display(['GA : fraction of positive gain = ' num2str(positive_gain_GA*100,'%.1f %%')]);
display(['GA gain = ' num2str(mean(gain_GA)*100,'%.1f %%')]);

positive_gain_GA_online = zeros(n_draw_done,1);
positive_gain_GA_online(gain_GA_online>=0) = 1/n_draw_done;
positive_gain_GA_online = sum(positive_gain_GA_online);

display(['GA online : fraction of positive gain = ' num2str(positive_gain_GA_online*100,'%.1f %%')]);
display(['GA online gain = ' num2str(mean(gain_GA_online)*100,'%.1f %%')]);

save(sprintf('online_computations_matlab_%d_particles_%d_draws',n_part,n_draw),'mse_reg','mse_GA','mse_GA_online','gain_GA','gain_GA_online')
%% 3. Display histogram 
if false 
figure;    
histogram(gain_GA,200,'Normalization','pdf')
xlim([-2, 1]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA algorithm')    
    
figure;    
histogram(gain_GA_online,200,'Normalization','pdf')
xlim([-2, 1]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA online algorithm')
end