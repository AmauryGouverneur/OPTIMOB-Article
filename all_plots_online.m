%% Plot function 
close all;

n_measurements = 10; %number of observations 
T = 25; %number of time steps 
measurements_spacing = 1; 
n_part_fine = 1000; 


n_part = 150; %number of particles 250
n_draw = 100; %number of draws for the MC MSE estimator 10

pop_size = 20; %population size of the GA algorithm 
max_gen = 20; %maximum number of generations

n_draw_comp = 500; %number of draws to realize 

try
    load(sprintf('here_new_GA_online_computations_%d_particles_%d_draws',n_part,n_draw));
    n_draw_done = length(find(mse_reg));
    mse_reg = [mse_reg(1:n_draw_done),zeros(1,n_draw_comp)];
    mse_GA = [mse_GA(1:n_draw_done),zeros(1,n_draw_comp)];
    mse_GA_online = [mse_GA_online(1:n_draw_done),zeros(1,n_draw_comp)];

catch
    mse_reg =  zeros(1,n_draw_comp); 
    mse_GA =  zeros(1,n_draw_comp); 
    mse_GA_online = zeros(1,n_draw_comp); 
    n_draw_done = 0;
end


%% 1. Comparison of the performances of the methods by running n_draw_comp draws
display(['Comparison of the performances of the methods by running draw number ',num2str(n_draw_done+1,'%.0f'),' to draw number ',num2str(n_draw_done + n_draw_comp,'%.0f')]);
t_start_comparison = tic;

meas_reg = round(linspace(0,T,n_measurements));

measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;


for j = n_draw_done+1:n_draw_done+n_draw_comp
x0 = initialization(1,0,0);
x = model(T,x0,0); %True state vector we want to reconstruct using online OIMPF
y = measurements(x,T,0);
% Computation of the measurement times reg, GF, RT 

disp('GA online (+GF) computations started');
t_start_GA_online = tic;
[meas_GA_online,meas_GA] = online_optimob(y,n_measurements,T,pop_size,max_gen,n_part,n_draw,measurements_spacing);
t_elapsed_GA_online = toc(t_start_GA_online);

measurements_GA = zeros(1,T+1); 
measurements_GA(meas_GA+1) = 1;

measurements_GA_online = zeros(1,T+1); 
measurements_GA_online(meas_GA_online+1) = 1;

part = initialization(n_part_fine,0,0);
tau_j_reg = particle_filter(y,measurements_reg,T,part,0);
tau_j_GA = particle_filter(y,measurements_GA,T,part,0);
tau_j_GA_online = particle_filter(y,measurements_GA_online,T,part,0);

mse_reg(j) = mean((objective(x,0)-tau_j_reg).^2,2);
mse_GA(j) =  mean((objective(x,0)-tau_j_GA).^2,2);
mse_GA_online(j) =  mean((objective(x,0)-tau_j_GA_online ).^2,2);

display(['GA online computations completed, mean gain =  ',num2str( (mean(mse_reg)-mean(mse_GA_online))/mean(mse_reg)*100,'%.1f %%'), 'gain apriori = ',num2str((mean(mse_reg)-mean(mse_GA))/mean(mse_reg)*100,'%.1f %%')  ,'time elapsed = ',num2str(t_elapsed_GA_online,'%.0f sec')]);

save(sprintf('here_new_GA_online_computations_%d_particles_%d_draws',n_part,n_draw),'mse_reg','mse_GA','mse_GA_online')
end
%% 2. Results and save

n_draw_done = length(find(mse_reg));
mse_reg = mse_reg(1:n_draw_done);
mse_GA = mse_GA(1:n_draw_done);
mse_GA_online = mse_GA_online(1:n_draw_done);

gain_GA = (mse_reg-mse_GA)./mse_reg;
gain_GA_online  = (mse_reg-mse_GA_online)./mse_reg; 



t_elapsed_comparison = toc(t_start_comparison);

display(['Computations of ',num2str(n_draw_done, '%.0f'),'draws, time elapsed = ',num2str(t_elapsed_comparison,'%.0f sec')]);

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

%% 3. Display histogram 
if true 
figure;    
histogram(gain_GA,50,'Normalization','pdf')
xlim([-2, 1]);
ylim([0,1.5]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA algorithm')    
    
figure;    
histogram(gain_GA_online,50,'Normalization','pdf')
xlim([-2, 1]);
ylim([0,1.5]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA online algorithm')
end