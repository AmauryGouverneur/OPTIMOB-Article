%% Plot function 
close all;
clc;
fontsize = 7*2;

n_draw_comp = 100000; %number of draws to build the histogram
n_part_fine = 1000;
n_measurements = 21; %number of observations 
T = 60; %number of time steps 

n_part = 200; %number of particles 250
n_draw = 100; %number of draws for the MC MSE estimator 100

pop_size = 45; %population size of the GA algorithm 
pop_size_SA = pop_size;
pop_size_GA = round(pop_size/2)*2;
max_gen = 25;
n_eval = pop_size*max_gen;

n_draw_GA = n_draw;
n_draw_SA = round(n_draw/2);

%GB_eval_max = (T+1);
%GB_eval_min = (n_measurements+1);
%n_eval_GB = (GB_eval_max+GB_eval_min)*(GB_eval_max-GB_eval_min+1)/2;
%n_average_meas_budget_GB = (T + n_measurements)/2;
%ratio_GB_draw = n_eval/n_eval_GB;
%ratio_GB_meas_budget = n_measurements/n_average_meas_budget_GB;
%n_part_GB = n_part;
%n_draw_GB = round(n_draw*ratio_GB_draw*ratio_GB_meas_budget);

GF_eval_max = (T+1);
GF_eval_min = (GF_eval_max-n_measurements+1);
n_eval_GF = (GF_eval_max+GF_eval_min)*(GF_eval_max-GF_eval_min+1)/2;
n_average_meas_budget_GF = (1 + n_measurements)/2;
ratio_GF_draw = n_eval/n_eval_GF; 
ratio_GF_meas_budget = n_measurements/n_average_meas_budget_GF;
n_part_GF = n_part;
n_draw_GF = round(n_draw*ratio_GF_draw*sqrt(ratio_GF_meas_budget));

GF_plot_spacing = linspace(0.1,1,10);
%% 1. Computation of the measurement times reg, GA, RT, SA, GF, GB 
display('1. Computation of the measurement times reg, GA, RT');

meas_reg = round(linspace(0,T,n_measurements));

%display('GB computations started'); 
%t_start_GB = tic;
%[meas_GB,~,avgCostHist_GB ,minCostHist_GB] = greedy_backward_algo(n_measurements,T,n_part_GB,n_draw_GB);
%t_elapsed_GB = toc(t_start_GB);
%display(['GB computations completed, time elapsed = ',num2str(t_elapsed_GB,'%.0f sec')]);
meas_GF = zeros(length(GF_plot_spacing),n_measurements);
avgCostHist_GF = zeros(length(GF_plot_spacing),n_measurements);
minCostHist_GF = zeros(length(GF_plot_spacing),n_measurements);

for i = 1:length(GF_plot_spacing)
    display('GF computations started');
    t_start_GF = tic;
    [meas_GF(i,:),~,avgCostHist_GF(i,:) ,minCostHist_GF(i,:)] = greedy_forward_algo(n_measurements,T,n_part_GF,round(n_draw_GF*GF_plot_spacing(i)));
    t_elapsed_GF = toc(t_start_GF);
    display(['GF computations completed, time elapsed = ',num2str(t_elapsed_GF,'%.0f sec')]);
end

display('SA computations started');
t_start_SA = tic;
[meas_SA,~,avgCostHist_SA ,minCostHist_SA] = SA_algo(n_measurements,T,pop_size_SA,max_gen,n_part,n_draw_SA);
t_elapsed_SA = toc(t_start_SA);
display(['SA computations completed, time elapsed = ',num2str(t_elapsed_SA,'%.0f sec')]);

display('GA computations started');
t_start_GA = tic;
[meas_GA,~,avgCostHist_GA ,minCostHist_GA] = genetical_algo(n_measurements,T,pop_size_GA,max_gen,n_part,n_draw_GA);
t_elapsed_GA = toc(t_start_GA);
display(['GA computations completed, time elapsed = ',num2str(t_elapsed_GA,'%.0f sec')]);

display('RT computations started');
t_start_RT = tic;
[meas_RT,~,avgCostHist_RT,minCostHist_RT] = random_trials(T,n_measurements, n_eval, n_part, n_draw);
t_elapsed_RT = toc(t_start_RT);
display(['RT computations completed, time elapsed = ',num2str(t_elapsed_RT,'%.0f sec')]);


%% 2. Comparison of the performances of the methods by running n_draw_comp draws
display(['2. Comparison of the performances of the methods by running ',num2str(n_draw_comp,'%.0f draws')]);
display('Comparison computations started');

t_start_comparison = tic;

measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;

measurements_GF = zeros(1,T+1); 
measurements_GF(meas_GF+1) = 1;

measurements_GB = zeros(1,T+1); 
measurements_GB(meas_GB+1) = 1;

measurements_SA = zeros(1,T+1); 
measurements_SA(meas_SA+1) = 1;

measurements_GA = zeros(1,T+1); 
measurements_GA(meas_GA+1) = 1;

measurements_RT = zeros(1,T+1); 
measurements_RT(meas_RT+1) = 1;

%0. Constants definition 

mse_reg = 0;
mse_GF = 0; 
mse_GB = 0; 
mse_SA = 0;
mse_GA = 0; 
mse_RT = 0;

gains_GF = zeros(n_draw_comp,1);
gains_GB = zeros(n_draw_comp,1);
gains_SA = zeros(n_draw_comp,1);
gains_GA = zeros(n_draw_comp,1);
gains_RT = zeros(n_draw_comp,1);

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
        tau_j_GB = particle_filter(y_j,measurements_GB,T,part,0);
        tau_j_SA = particle_filter(y_j,measurements_SA,T,part,0);
        tau_j_GA = particle_filter(y_j,measurements_GA,T,part,0);
        tau_j_RT = particle_filter(y_j,measurements_RT,T,part,0);
        
        %3. MSE computation
        err_reg = mean((objective(x_j)-objective(tau_j_reg)).^2,2);
        err_GF =  mean((objective(x_j)-objective(tau_j_GF)).^2,2);
        err_GB = mean((objective(x_j)-objective(tau_j_GB)).^2,2);
        err_SA =  mean((objective(x_j)-objective(tau_j_SA)).^2,2);
        err_GA =  mean((objective(x_j)-objective(tau_j_GA)).^2,2);
        err_RT = mean((objective(x_j)-objective(tau_j_RT)).^2,2);
        
        mse_reg = mse_reg + 1/n_draw_comp*err_reg;
        mse_GF = mse_GF + 1/n_draw_comp*err_GF; 
        mse_GB = mse_GB + 1/n_draw_comp*err_GB;
        mse_SA = mse_SA + 1/n_draw_comp*err_SA; 
        mse_GA = mse_GA + 1/n_draw_comp*err_GA; 
        mse_RT = mse_RT + 1/n_draw_comp*err_RT;
        
        gains_GF(j) = (err_reg - err_GF)/err_reg;
        gains_GB(j) = (err_reg - err_GB)/err_reg;
        gains_SA(j) = (err_reg - err_SA)/err_reg;
        gains_GA(j) = (err_reg - err_GA)/err_reg;
        gains_RT(j) = (err_reg - err_RT)/err_reg;
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
histogram(gains_GA,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2.2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GF algorithm')

positive_gain_GB = zeros(n_draw_comp,1);
positive_gain_GB(gains_GB>=0) = 1/n_draw_comp;
positive_gain_GB = sum(positive_gain_GB);
display(['GB : average gain = ' num2str(mean(gains_GB)*100,'%.1f %%')]);
display(['GB : fraction of positive gain = ' num2str(positive_gain_GB*100,'%.1f %%')]);
figure;    
histogram(gains_GB,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2.2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GB algorithm')

positive_gain_SA = zeros(n_draw_comp,1);
positive_gain_SA(gains_SA>=0) = 1/n_draw_comp;
positive_gain_SA = sum(positive_gain_SA);
display(['SA : average gain = ' num2str(mean(gains_SA)*100,'%.1f %%')]);
display(['SA : fraction of positive gain = ' num2str(positive_gain_SA*100,'%.1f %%')]);
figure;    
histogram(gains_SA,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2.2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with SA algorithm')

positive_gain_GA = zeros(n_draw_comp,1);
positive_gain_GA(gains_GA>=0) = 1/n_draw_comp;
positive_gain_GA = sum(positive_gain_GA);
display(['GA : average gain = ' num2str(mean(gains_GA)*100,'%.1f %%')]);
display(['GA : fraction of positive gain = ' num2str(positive_gain_GA*100,'%.1f %%')]);
figure;    
histogram(gains_GA,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2.2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA algorithm')

positive_gain_RT = zeros(n_draw_comp,1);
positive_gain_RT(gains_RT>=0) = 1/n_draw_comp;
positive_gain_RT = sum(positive_gain_RT);
display(['RT : average gain = ' num2str(mean(gains_RT)*100,'%.1f %%')]);
display(['RT : fraction of positive gain = ' num2str(positive_gain_RT*100,'%.1f %%')]);
figure;
histogram(gains_RT,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2.2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with RT algorithm')

%%
mse_reg = 31.076762302039494; 
FIG = figure; hold on;
genList = (1:length(avgCostHist_GA))-1;
evalList_RT = (1:n_eval);
evalList_GF = cumsum((GF_eval_max:-1:GF_eval_min)*n_eval_GF);
evalList_GF = round(evalList_GF*(n_eval/evalList_GF(end)));
evalList_GB = cumsum((GB_eval_max:-1:GB_eval_min)*n_eval_GB);
evalList_GB = round(evalList_GB*(n_eval/evalList_GB(end)));

plot(round(genList*pop_size)*n_draw,0*avgCostHist_GA+mse_reg,'-.k');
plot(round(genList*pop_size)*n_draw,avgCostHist_GA,'--r','LineWidth',1.5);
plot(round(genList*pop_size)*n_draw,minCostHist_GA,'-r','LineWidth',1.5);
plot(round(genList*pop_size)*n_draw,avgCostHist_SA,'--m');
plot(round(genList*pop_size)*n_draw,minCostHist_SA,'-m');
plot(evalList_RT*n_draw,avgCostHist_RT,'--k');
plot(evalList_RT*n_draw,minCostHist_RT,'-k');


plot(round(evalList_GF(end)*n_draw*GF_plot_spacing),avgCostHist_GF(:,end),'--b');
plot(round(evalList_GF(end)*n_draw*GF_plot_spacing),minCostHist_GF(:,end),'-b');


%plot(evalList_GB,avgCostHist_GB,'--c');
%plot(evalList_GB,minCostHist_GB,'-c');
box on
xlim([0 n_eval*n_draw])
ylim([20 40])
legend('regualar MSE cost','GA: average cost', 'GA: min cost GA','SA: average cost', 'SA: min cost GA','RT: average cost','RT: min cost','GF: average cost','GF: min cost')
xlabel('number of draws','interpreter','latex')
ylabel('$\hat{\mathrm{E}}_{\mathrm{MSE}}$','interpreter','latex')
title('Evolution of the average and minimum cost $\hat{\mathrm{E}}_{\mathrm{MSE}}$','interpreter','latex')
set(findall(gcf,'-property','FontSize'),'Fontsize',fontsize);

