%% Plot function 
close all;
fontsize = 7*2;
%Rectify the computations times to have something to compare
n_draw_comp = 1000; %number of draws to build the histogram
n_part_fine = 200;
n_measurements =10; %number of observations 
T = 25; %number of time steps 
measurements_spacing = 1; 

to_test = [1,0,1,1,0]; %GF,GB,SA,GA,RT

n_part = 150; %number of particles 250
n_draw = 100; %number of draws for the MC MSE estimator 100

pop_size = 20; %population size of the GA algorithm 
pop_size_SA = round(pop_size/2);
pop_size_GA = 2*round(pop_size/2);
max_gen = 20;
max_gen_SA = max_gen;

n_eval = pop_size*max_gen;

n_draw_GA = n_draw;
n_draw_SA = n_draw;

GB_eval_max = (T+1);
GB_eval_min = (n_measurements+1);
n_eval_GB = (GB_eval_max+GB_eval_min)*(GB_eval_max-GB_eval_min+1)/2;
n_average_meas_budget_GB = (T + n_measurements)/2;
ratio_GB_draw = n_eval/n_eval_GB;
ratio_GB_meas_budget = n_measurements/n_average_meas_budget_GB;
n_part_GB = n_part;
n_draw_GB = round(n_draw*ratio_GB_draw);

GF_eval_max = (T+1);
GF_eval_min = (GF_eval_max-n_measurements+1);
n_eval_GF = (GF_eval_max+GF_eval_min)*(GF_eval_max-GF_eval_min+1)/2;
n_average_meas_budget_GF = (1 + n_measurements)/2;
ratio_GF_draw = n_eval/n_eval_GF; 
ratio_GF_meas_budget = n_measurements/n_average_meas_budget_GF;
n_part_GF = n_part;
n_draw_GF = round(n_draw*ratio_GF_draw);
%n_draw_GF = 150;

%Greedy_plot_spacing = linspace(0.2,1,5);
Greedy_plot_spacing = [1];

mean_gain_GF = 0; 
mean_gain_GB = 0;
mean_gain_SA = 0;
mean_gain_GA = 0; 
mean_gain_RT = 0;

n_test = 10;
for l = 1:n_test
%% 1. Computation of the measurement times reg, GA, RT, SA, GF, GB 
disp('1. Computation of the measurement times reg, GA, RT, SA, GF, GB');

meas_reg = round(linspace(0,T,n_measurements));


meas_GF = zeros(length(Greedy_plot_spacing),n_measurements);
minCostEnd_GF = zeros(1,length(Greedy_plot_spacing));

if to_test(1)
for i = 1:length(Greedy_plot_spacing)
    disp('GF computations started');
    t_start_GF = tic;
    [meas_GF(i,:),minCostEnd_GF(i),~,~] = greedy_forward_algo(n_measurements,T,n_part_GF,round(n_draw_GF*Greedy_plot_spacing(i)),measurements_spacing);
    t_elapsed_GF = toc(t_start_GF);
    display(['GF computations completed, time elapsed = ',num2str(t_elapsed_GF,'%.0f sec')]);
    
end
%%
end
if to_test(2)
meas_GB = zeros(length(Greedy_plot_spacing),n_measurements);
minCostEnd_GB = zeros(1,length(Greedy_plot_spacing));
for i = 1:length(Greedy_plot_spacing)   
    disp('GB computations started'); 
    t_start_GB = tic;
    [meas_GB(i,:),minCostEnd_GB(i),~,~] = greedy_backward_algo(n_measurements,T,n_part_GB,round(n_draw_GB*Greedy_plot_spacing(i)),measurements_spacing);
    t_elapsed_GB = toc(t_start_GB);
    display(['GB computations completed, time elapsed = ',num2str(t_elapsed_GB,'%.0f sec')]);
end
end
%%
if to_test(3)
disp('SA computations started');
t_start_SA = tic;
[meas_SA,~,avgCostHist_SA ,minCostHist_SA] = SA_algo(n_measurements,T,pop_size_SA,max_gen_SA,n_part,n_draw_SA);
t_elapsed_SA = toc(t_start_SA);
display(['SA computations completed, time elapsed = ',num2str(t_elapsed_SA,'%.0f sec')]);
%%
end
if to_test(4)
disp('GA computations started');
t_start_GA = tic;
[meas_GA,~,avgCostHist_GA ,minCostHist_GA] = genetical_algo(n_measurements,T,pop_size_GA,max_gen,n_part,n_draw_GA,measurements_spacing);
t_elapsed_GA = toc(t_start_GA);
display(['GA computations completed, time elapsed = ',num2str(t_elapsed_GA,'%.0f sec')]);
%%
end
if to_test(5)
disp('RT computations started');
t_start_RT = tic;
[meas_RT,~,avgCostHist_RT,minCostHist_RT] = random_trials(T,n_measurements, n_eval, n_part, n_draw,measurements_spacing);
t_elapsed_RT = toc(t_start_RT);
display(['RT computations completed, time elapsed = ',num2str(t_elapsed_RT,'%.0f sec')]);
end

%% 2. Comparison of the performances of the methods by running n_draw_comp draws
display(['2. Comparison of the performances of the methods by running ',num2str(n_draw_comp,'%.0f draws')]);
disp('Comparison computations started');

t_start_comparison = tic;

measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;
mse_reg = 0;


if to_test(1)
    measurements_GF = zeros(1,T+1); 
    measurements_GF(meas_GF(end,:)+1) = 1;
    mse_GF = 0; 
    gains_GF = zeros(n_draw_comp,1);
end
if to_test(2)
    measurements_GB = zeros(1,T+1); 
    measurements_GB(meas_GB(end,:)+1) = 1;
    mse_GB = 0; 
    gains_GB = zeros(n_draw_comp,1);
end
if to_test(3)
    measurements_SA = zeros(1,T+1); 
    measurements_SA(meas_SA+1) = 1;
    mse_SA = 0;
    gains_SA = zeros(n_draw_comp,1);
end
if to_test(4)
    measurements_GA = zeros(1,T+1); 
    measurements_GA(meas_GA+1) = 1;
    mse_GA = zeros(n_draw_comp,1);
    gains_GA = zeros(n_draw_comp,1);
end
if to_test(5)
    measurements_RT = zeros(1,T+1); 
    measurements_RT(meas_RT+1) = 1;
    mse_RT = 0;
    gains_RT = zeros(n_draw_comp,1);
end
%0. Constants definition 



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
        err_reg = mean((objective(x_j,0)-tau_j_reg).^2,2);
        mse_reg = mse_reg + 1/n_draw_comp*err_reg;
        
        if to_test(1)
            tau_j_GF = particle_filter(y_j,measurements_GF,T,part,0);
            err_GF =  mean((objective(x_j,0)-tau_j_GF).^2,2);
            mse_GF = mse_GF + 1/n_draw_comp*err_GF;
            gains_GF(j) = (err_reg - err_GF)/err_reg;
        end
        if to_test(2)
            tau_j_GB = particle_filter(y_j,measurements_GB,T,part,0);
            err_GB = mean((objective(x_j,0)-tau_j_GB).^2,2);
            mse_GB = mse_GB + 1/n_draw_comp*err_GB;
            gains_GB(j) = (err_reg - err_GB)/err_reg;
        end
        if to_test(3)
            tau_j_SA = particle_filter(y_j,measurements_SA,T,part,0);
            err_SA =  mean((objective(x_j,0)-tau_j_SA).^2,2);
            mse_SA = mse_SA + 1/n_draw_comp*err_SA;
            gains_SA(j) = (err_reg - err_SA)/err_reg;
        end
        if to_test(4)
            tau_j_GA = particle_filter(y_j,measurements_GA,T,part,0);
            err_GA =  mean((objective(x_j,0)-tau_j_GA).^2,2);
            mse_GA(j) = err_GA;
            gains_GA(j) = (err_reg - err_GA)/err_reg;
        end
        if to_test(5)
            tau_j_RT = particle_filter(y_j,measurements_RT,T,part,0);
            err_RT = mean((objective(x_j,0)-tau_j_RT).^2,2); 
            mse_RT = mse_RT + 1/n_draw_comp*err_RT;
            gains_RT(j) = (err_reg - err_RT)/err_reg;
        end
             
end

t_elapsed_comparison = toc(t_start_comparison);
display(['Comparison computations completed, time elapsed = ',num2str(t_elapsed_comparison,'%.0f sec')]);
%%
if to_test(1)
positive_gain_GF = zeros(n_draw_comp,1);
positive_gain_GF(gains_GF>=0) = 1/n_draw_comp;
positive_gain_GF = sum(positive_gain_GF);
display(['GF : average gain = ' num2str(mean(gains_GF)*100,'%.1f %%') ' (fraction of positive gain = ' num2str(positive_gain_GF*100,'%.1f %%)')]);
mean_gain_GF = mean_gain_GF + mean(gains_GF)*100/n_test ; 
%display();
if 0 
figure;    
histogram(gains_GF,200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GF algorithm')
end
end 
if to_test(2)
positive_gain_GB = zeros(n_draw_comp,1);
positive_gain_GB(gains_GB>=0) = 1/n_draw_comp;
positive_gain_GB = sum(positive_gain_GB);
display(['GB : average gain = ' num2str(mean(gains_GB)*100,'%.1f %%') ' (fraction of positive gain = ' num2str(positive_gain_GB*100,'%.1f %%)')]);
mean_gain_GB = mean_gain_GB + mean(gains_GB)*100/n_test ; 
if 0
figure;    
histogram(gains_GB,200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GB algorithm')
end
end
if to_test(3)
positive_gain_SA = zeros(n_draw_comp,1);
positive_gain_SA(gains_SA>=0) = 1/n_draw_comp;
positive_gain_SA = sum(positive_gain_SA);
display(['SA : average gain = ' num2str(mean(gains_SA)*100,'%.1f %%') ' (fraction of positive gain = ' num2str(positive_gain_SA*100,'%.1f %%)')]);
mean_gain_SA = mean_gain_SA + mean(gains_SA)*100/n_test ; 
if 0
figure;    
histogram(gains_SA,200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with SA algorithm')
end
end
if to_test(4)
positive_gain_GA = zeros(n_draw_comp,1);
positive_gain_GA(gains_GA>=0) = 1/n_draw_comp;
positive_gain_GA = sum(positive_gain_GA);
display(['GA : average gain = ' num2str(mean(gains_GA)*100,'%.1f %%') ' (fraction of positive gain = ' num2str(positive_gain_GA*100,'%.1f %%)')]);
mean_gain_GA = mean_gain_GA + mean(gains_GA)*100/n_test ; 
if 0
figure;    
histogram(gains_GA,200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA algorithm')
end
end
if to_test(5)
positive_gain_RT = zeros(n_draw_comp,1);
positive_gain_RT(gains_RT>=0) = 1/n_draw_comp;
positive_gain_RT = sum(positive_gain_RT);
display(['RT : average gain = ' num2str(mean(gains_RT)*100,'%.1f %%') ' (fraction of positive gain = ' num2str(positive_gain_RT*100,'%.1f %%)')]);
mean_gain_RT = mean_gain_RT + mean(gains_RT)*100/n_test ; 
if 0
figure;
histogram(gains_RT,200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with RT algorithm')
end
end
end



%%
%mse_reg = 31.076762302039494; 
if sum(to_test/5)==1
FIG = figure; hold on;
genList = (1:length(avgCostHist_GA))-1;
evalList_RT = (1:n_eval);
evalList_GF = cumsum((GF_eval_max:-1:GF_eval_min)*n_eval_GF);
evalList_GF = round(evalList_GF*(n_eval/evalList_GF(end)));
evalList_GB = cumsum((GB_eval_max:-1:GB_eval_min)*n_eval_GB);
evalList_GB = round(evalList_GB*(n_eval/evalList_GB(end)));
the_end = round(genList*pop_size) ;
plot(round(genList*pop_size),0*avgCostHist_GA+mse_reg,'-.k');
plot(round(genList*pop_size),avgCostHist_GA,'--r','LineWidth',1.5);
plot(round(genList*pop_size),minCostHist_GA,'-r','LineWidth',1.5);
%plot(round(genList*pop_size)*n_draw,avgCostHist_SA,'--m');
%plot((0:length(minCostHist_SA)-1)*the_end(end)/(length(minCostHist_SA)-1),minCostHist_SA,'-m');
plot(round(genList*pop_size),minCostHist_SA,'-m');

%plot(evalList_RT*n_draw,avgCostHist_RT,'--k');
plot(evalList_RT,minCostHist_RT,'-k');


plot(round(evalList_GF(end)*Greedy_plot_spacing),minCostEnd_GF(:),'*b','LineWidth',1.5);
plot(round(evalList_GB(end)*Greedy_plot_spacing),minCostEnd_GB(:),'*c','LineWidth',1.5);

box on
xlim([0 n_eval])
ylim([10*floor(minCostHist_GA(end)/10) 10*ceil(minCostHist_RT(1)/10)])
%legend('regualar MSE cost','GA: average cost', 'GA: min cost GA', 'SA: min cost GA','RT: min cost','GF: min cost')
legend('regualar MSE cost','GA: average cost', 'GA: min cost GA', 'SA: min cost SA','RT: min cost','GF: min cost','GB: min cost')
%legend('regualar MSE cost','GA: average cost', 'GA: min cost GA','SA: average cost', 'SA: min cost GA','RT: average cost','RT: min cost','GF: average cost','GF: min cost')
xlabel('number of cost function evaluation','interpreter','latex')
ylabel('$\hat{\mathrm{E}}_{\mathrm{MSE}}$','interpreter','latex')
title('Evolution of the average and minimum cost $\hat{\mathrm{E}}_{\mathrm{MSE}}$','interpreter','latex')
set(findall(gcf,'-property','FontSize'),'Fontsize',fontsize);
end

%%
if to_test(4)
disp('Comparison over a single draw : GF and ref ');
meas_reg = round(linspace(0,T,n_measurements));

measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;

measurements_GA = zeros(1,T+1); 
measurements_GA(meas_GA+1) = 1;
figure
x0 = initialization(1);
x_j = model(T,x0);
y_j = measurements(x_j,T,0);
part = initialization(n_part_fine,0,0);
tau_j_reg = particle_filter(y_j,measurements_reg,T,part,0);
tau_j_GA = particle_filter(y_j,measurements_GA,T,part,0);
% err_reg = mean((objective(x_j,0)-objective(tau_j_reg,0)).^2,2);
% err_GF = mean((objective(x_j,0)-objective(tau_j_GF,0)).^2,2);
err_reg = mean((objective(x_j,0)-tau_j_reg).^2,2);
err_GA = mean((objective(x_j,0)-tau_j_GA).^2,2);

gain_GA = (err_reg - err_GA)/err_reg;
display(['GA : gain = ' num2str((gain_GA)*100,'%.1f %%')]);

plot(0:T,(objective(x_j,0)),'k')
hold on
% plot(0:T,objective(tau_j_reg,0),':b')
% plot(0:T,objective(tau_j_GF,0),'-.r')
plot(0:T,tau_j_reg,':b')
plot(0:T,tau_j_GA,'-.r')
y_limits_plot = ylim;
plot(meas_reg,y_limits_plot(2)*ones(n_measurements,1),'+b')
plot(meas_GA,y_limits_plot(1)*ones(n_measurements,1),'*r')
legend('truth','reg','optimized')
ylabel('z(t)','interpreter','latex')
xlabel('time, t','interpreter','latex')
xlim([0 T])

end
%%
if 0 
disp('Samples');

a_minus = 8.8;
a_plus = 24;

omega_minus = 0.13;
omega_plus = 0.21;

b_minus = -5.8;
b_plus = 5.8;
x0 = initialization(1);
x_j = model(T,x0);

figure
hold on

subplot(4,1,[1,2]);
plot(linspace(0,T/5,T+1),(objective(x_j,0)),'k')
ylabel('motion [mm] ','interpreter','latex')

subplot(4,1,3);
hold on
ylabel('amplitude [mm] ','interpreter','latex')
plot(linspace(0,T/5,T+1),x_j(1,:),'k')
plot(linspace(0,T/5,T+1),a_minus*ones(1,T+1))
plot(linspace(0,T/5,T+1),a_plus*ones(1,T+1))

subplot(4,1,4);
hold on
ylabel('offset [mm] ','interpreter','latex')
plot(linspace(0,T/5,T+1),b_minus*ones(1,T+1))
plot(linspace(0,T/5,T+1),b_plus*ones(1,T+1))
plot(linspace(0,T/5,T+1),x_j(3,:),'k')


xlabel('time, [s]','interpreter','latex')
set(findall(gcf,'-property','FontSize'),'Fontsize',fontsize);
end