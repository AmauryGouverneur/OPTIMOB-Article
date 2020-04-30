fontsize = 7*2;

pop_size = 50; %population size of the GA algorithm 
max_gen = 20;

n_eval = pop_size*max_gen;

n_draw_comp = 10000; %number of draws to build the histogram
n_part_fine = 1000;

n_measurements =11; %number of observations 
T = 30; %number of time steps 

n_part = 200; %number of particles 250
n_draw = 200; %number of draws for the MC MSE estimator 100

GFR_eval_max = (T+1);
GFR_eval_min = (GFR_eval_max-n_measurements+1);
n_eval_GFR = (GFR_eval_max+GFR_eval_min)*(GFR_eval_max-GFR_eval_min+1)/2+GFR_eval_min*n_measurements;
n_average_meas_budget_GFR = 0.5*((1 + n_measurements)/2+n_measurements);
ratio_GFR_draw = n_eval/n_eval_GFR; 
ratio_GFR_meas_budget = n_measurements/n_average_meas_budget_GFR;
n_part_GFR = n_part;
n_draw_GFR = round(n_draw*ratio_GFR_draw*sqrt(ratio_GFR_meas_budget))

GF_eval_max = (T+1);
GF_eval_min = (GF_eval_max-n_measurements+1);
n_eval_GF = (GF_eval_max+GF_eval_min)*(GF_eval_max-GF_eval_min+1)/2;
n_average_meas_budget_GF = (1 + n_measurements)/2;
ratio_GF_draw = n_eval/n_eval_GF; 
ratio_GF_meas_budget = n_measurements/n_average_meas_budget_GF;
n_part_GF = n_part;
n_draw_GF = round(n_draw*ratio_GF_draw*sqrt(ratio_GF_meas_budget))
%%
display('GF computations started');
t_start_GF = tic;
[meas_GF,minCostEnd_GF,~,~] = greedy_forward_algo(n_measurements,T,n_part_GF,n_draw_GF);
t_elapsed_GF = toc(t_start_GF);
display(['GF computations completed, time elapsed = ',num2str(t_elapsed_GF,'%.0f sec')]);
%%
display('GF release computations started');
t_start_GFR = tic;
[meas_GFR,minCostEnd_GFR,~,~] = greedy_forward_with_release_algo(n_measurements,T,n_part_GFR,n_draw_GFR);
t_elapsed_GFR = toc(t_start_GFR);

%%
display(['GFR computations completed, time elapsed = ',num2str(t_elapsed_GFR,'%.0f sec')]);

display(['2. Comparison of the performances of the methods by running ',num2str(n_draw_comp,'%.0f draws')]);
display('Comparison computations started');

t_start_comparison = tic;

meas_reg = round(linspace(0,T,n_measurements));


measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;

measurements_GF = zeros(1,T+1); 
measurements_GF(meas_GF(end,:)+1) = 1;

measurements_GFR = zeros(1,T+1); 
measurements_GFR(meas_GFR(end,:)+1) = 1;

%0. Constants definition 

mse_reg = 0;
mse_GF = 0; 
mse_GFR = 0; 


gains_GF = zeros(n_draw_comp,1);
gains_GFR = zeros(n_draw_comp,1);


x0 = initialization(n_draw_comp,0,0);
for j = 1:n_draw_comp
        %1.Simulation

        %1.1. Random motion model : X_j
        x_j = model(T,x0(j),0);

        %1.2. Artificial data record Y_j
        y_j = measurements(x_j,T,0);

        %2. Filtering
        part = initialization(n_part_fine,0,0);
        tau_j_reg = particle_filter(y_j,measurements_reg,T,part,0);
        tau_j_GF = particle_filter(y_j,measurements_GF,T,part,0);
        tau_j_GFR = particle_filter(y_j,measurements_GFR,T,part,0);

        
        %3. MSE computation
        err_reg = mean((objective(x_j,0)-objective(tau_j_reg,0)).^2,2);
        err_GF =  mean((objective(x_j,0)-objective(tau_j_GF,0)).^2,2);
        err_GFR =  mean((objective(x_j,0)-objective(tau_j_GFR,0)).^2,2);

       
        mse_reg = mse_reg + 1/n_draw_comp*err_reg;
        mse_GF = mse_GF + 1/n_draw_comp*err_GF; 
        mse_GFR = mse_GFR + 1/n_draw_comp*err_GFR; 

        gains_GF(j) = (err_reg - err_GF)/err_reg;
        gains_GFR(j) = (err_reg - err_GFR)/err_reg;
end

t_elapsed_comparison = toc(t_start_comparison);
display(['Comparison computations completed, time elapsed = ',num2str(t_elapsed_comparison,'%.0f sec')]);
%%
positive_gain_GF = zeros(n_draw_comp,1);
positive_gain_GF(gains_GF>=0) = 1/n_draw_comp;
positive_gain_GF = sum(positive_gain_GF);
display(['GF : average gain = ' num2str(mean(gains_GF)*100,'%.1f %%')]);
display(['GF : fraction of positive gain = ' num2str(positive_gain_GF*100,'%.1f %%')]);
figure;    
histogram(gains_GF,'Normalization','pdf')
% xlim([-2, 1]);
% ylim([0, 2.2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GF algorithm')


positive_gain_GFR = zeros(n_draw_comp,1);
positive_gain_GFR(gains_GFR>=0) = 1/n_draw_comp;
positive_gain_GFR = sum(positive_gain_GFR);
display(['GFR : average gain = ' num2str(mean(gains_GFR)*100,'%.1f %%')]);
display(['GFR : fraction of positive gain = ' num2str(positive_gain_GFR*100,'%.1f %%')]);
figure;    
histogram(gains_GFR,'Normalization','pdf')
% xlim([-2, 1]);
% ylim([0, 2.2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GFR algorithm')

