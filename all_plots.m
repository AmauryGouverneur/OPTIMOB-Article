%% Plot function 
close all;
fontsize = 7*2.5;
%Rectify the computations times to have something to compare
n_draw_comp = 100000; %number of draws to build the histogram
n_part_fine = 1000;
n_measurements =11; %number of observations 
T = 30; %number of time steps 
measurements_spacing = 1; 

to_test = [0,0,1,1,0]; %GF,GB,SA,GA,RT

n_part = 200; %number of particles 250
n_draw = 1000; %number of draws for the MC MSE estimator 100

pop_size = 50; %population size of the GA algorithm 
pop_size_SA = round(pop_size/2);
pop_size_GA = 2*round(pop_size/2);
max_gen = 25;
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
%% 1. Computation of the measurement times reg, GA, RT, SA, GF, GB 


disp('1. Computation of the measurement times reg, GA, RT, SA, GF, GB');

meas_reg = round(linspace(0,T,n_measurements));




if to_test(1)
meas_GF = zeros(length(Greedy_plot_spacing),n_measurements);
minCostEnd_GF = zeros(1,length(Greedy_plot_spacing));
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
mse_reg = zeros(n_draw_comp,1);


if to_test(1)
    measurements_GF = zeros(1,T+1); 
    measurements_GF(meas_GF(end,:)+1) = 1;
    mse_GF = zeros(n_draw_comp,1);
    gains_GF = zeros(n_draw_comp,1);
end
if to_test(2)
    measurements_GB = zeros(1,T+1); 
    measurements_GB(meas_GB(end,:)+1) = 1;
    mse_GB = zeros(n_draw_comp,1); 
    gains_GB = zeros(n_draw_comp,1);
end
if to_test(3)
    measurements_SA = zeros(1,T+1); 
    measurements_SA(meas_SA+1) = 1;
    mse_SA = zeros(n_draw_comp,1);
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
    mse_RT = zeros(n_draw_comp,1);
    gains_RT = zeros(n_draw_comp,1);
end
%0. Constants definition 



x0 = initialization(n_draw_comp,0,0);
for j = 1:n_draw_comp
        %1.Simulation

        %1.1. Random motion model : X_j
        x_j = model(T,x0(:,j),0);

        %1.2. Artificial data record Y_j
        y_j = measurements(x_j,T,0);

        %2. Filtering
        part = initialization(n_part_fine,0,0);
        tau_j_reg = particle_filter(y_j,measurements_reg,T,part,0);
        err_reg = mean((objective(x_j,0)-tau_j_reg).^2);
        mse_reg(j) = err_reg;
        
        if to_test(1)
            tau_j_GF = particle_filter(y_j,measurements_GF,T,part,0);
            err_GF =  mean((objective(x_j,0)-tau_j_GF).^2);
            mse_GF(j) = err_GF;
            gains_GF(j) = (err_reg - err_GF)/err_reg;
        end
        if to_test(2)
            tau_j_GB = particle_filter(y_j,measurements_GB,T,part,0);
            err_GB = mean((objective(x_j,0)-tau_j_GB).^2);
            mse_GB(j) = err_GB;
            gains_GB(j) = (err_reg - err_GB)/err_reg;
        end
        if to_test(3)
            tau_j_SA = particle_filter(y_j,measurements_SA,T,part,0);
            err_SA =  mean((objective(x_j,0)-tau_j_SA).^2);
            mse_SA(j) = err_SA;
            gains_SA(j) = (err_reg - err_SA)/err_reg;
        end
        if to_test(4)
            tau_j_GA = particle_filter(y_j,measurements_GA,T,part,0);
            err_GA =  mean((objective(x_j,0)-tau_j_GA).^2);
            mse_GA(j) = err_GA;
            gains_GA(j) = (err_reg - err_GA)/err_reg;
        end
        if to_test(5)
            tau_j_RT = particle_filter(y_j,measurements_RT,T,part,0);
            err_RT = mean((objective(x_j,0)-tau_j_RT).^2); 
            mse_RT(j) = err_RT;
            gains_RT(j) = (err_reg - err_RT)/err_reg;
        end
             
end

t_elapsed_comparison = toc(t_start_comparison);
display(['Comparison computations completed, time elapsed = ',num2str(t_elapsed_comparison,'%.0f sec')]);
%%
alpha = 0.05; 
B = 1000; 

mse_error_margin_reg = norminv(1-alpha/2)*std(mse_reg)/(sqrt(n_draw_comp)*mean(mse_reg));
visualization_1 = 1;
visualization_2 = 0;

if to_test(1)
mean_gain_error_margin_GF = norminv(1-alpha/2)*std(gains_GF)/(sqrt(n_draw_comp));
mse_error_margin_GF = norminv(1-alpha/2)*std(mse_GF)/(sqrt(n_draw_comp)*mean(mse_GF));

bootstrap_error_median_GF = sort(median(gains_GF)-bootstrp(B,@median,gains_GF));
bootstrap_error_pos_frac_gain_GF = sort(mean(gains_GA>=0)-bootstrp(B,@(x) mean(x>=0),gains_GA));

median_error_margin_GF = [bootstrap_error_median_GF(ceil(alpha/2*B)),bootstrap_error_median_GF(ceil((1-alpha/2)*B))];
pos_frac_error_margin_GF = [bootstrap_error_pos_frac_gain_GF(ceil(alpha/2*B)),bootstrap_error_pos_frac_gain_GF(ceil((1-alpha/2)*B))];

display(['GF : gain on expected mse = ',num2str((1-(mean(mse_GF))/mean(mse_reg))*100,'%.1f%%'),' (+/- ',num2str((mse_error_margin_GF)/(1-mse_error_margin_reg)*100,'%.1f%%)') , ', average gain = ' num2str(mean(gains_GF)*100,'%.1f%%'),' (+/- ',num2str(mean_gain_error_margin_GF*100,'%.1f%%)'), ', median gain = ' num2str(mean((median(gains_GF)+median_error_margin_GF)*100),'%.1f%%'),' (+/- ',num2str(mean(abs(median_error_margin_GF))*100,'%.1f%%)'), ', fraction of positive gain = ' num2str(mean((mean(gains_GF>=0)+pos_frac_error_margin_GF)*100),'%.1f%%'),' (+/- ',num2str(mean(abs(pos_frac_error_margin_GF))*100,'%.1f%%)')]);

%display();
if visualization_1
figure;    
histogram(gains_GF,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 1.5]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GF algorithm')
end
if visualization_2
figure;    
histogram((mean(mse_reg)-mse_GF)/mean(mse_reg),200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('gain vs expected MSE$_{reg}$, $g$','interpreter','latex')
%title('Histogram of the gain vs expected MSE$_{reg}$ obtained with GF algorithm','interpreter','latex')
end
end 
if to_test(2)
mean_gain_error_margin_GB = norminv(1-alpha/2)*std(gains_GB)/(sqrt(n_draw_comp));
mse_error_margin_GB = norminv(1-alpha/2)*std(mse_GB)/(sqrt(n_draw_comp)*mean(mse_GB));

bootstrap_error_median_GB = sort(median(gains_GB )-bootstrp(B,@median,gains_GB ));
bootstrap_error_pos_frac_gain_GB = sort(mean(gains_GB >= 0)-bootstrp(B,@(x) mean(x>=0),gains_GB ));

median_error_margin_GB = [bootstrap_error_median_GB(ceil(alpha/2*B)),bootstrap_error_median_GB(ceil((1-alpha/2)*B))];
pos_frac_error_margin_GB = [bootstrap_error_pos_frac_gain_GB(ceil(alpha/2*B)),bootstrap_error_pos_frac_gain_GB(ceil((1-alpha/2)*B))];

display(['GB : gain on expected mse = ',num2str((1-(mean(mse_GB ))/mean(mse_reg))*100,'%.1f%%'),' (+/- ',num2str((mse_error_margin_GB )/(1-mse_error_margin_reg)*100,'%.1f%%)') , ', average gain = ' num2str(mean(gains_GB )*100,'%.1f%%'),' (+/- ',num2str(mean_gain_error_margin_GB*100,'%.1f%%)'), ', median gain = ' num2str(mean((median(gains_GB )+median_error_margin_GB )*100),'%.1f%%'),' (+/- ',num2str(mean(abs(median_error_margin_GB ))*100,'%.1f%%)'), ', fraction of positive gain = ' num2str(mean((mean(gains_GB>=0)+pos_frac_error_margin_GB )*100),'%.1f%%'),' (+/- ',num2str(mean(abs(pos_frac_error_margin_GB ))*100,'%.1f%%)')]);
if visualization_1
figure;    
histogram(gains_GB,200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 1.5]); 
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GB algorithm')
end
if visualization_2
figure;    
histogram((mean(mse_reg)-mse_GB)/mean(mse_reg),200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('gain vs expected MSE$_{reg}$, $g$','interpreter','latex')
title('Histogram of the gain vs expected MSE$_{reg}$ obtained with GB algorithm','interpreter','latex')
end
end
if to_test(3)
mean_gain_error_margin_SA = norminv(1-alpha/2)*std(gains_SA)/(sqrt(n_draw_comp));
mse_error_margin_SA = norminv(1-alpha/2)*std(mse_SA)/(sqrt(n_draw_comp)*mean(mse_SA));

bootstrap_error_median_SA = sort(median(gains_SA)-bootstrp(B,@median,gains_SA));
bootstrap_error_pos_frac_gain_SA = sort(mean(gains_SA >=0)-bootstrp(B,@(x) mean(x>=0),gains_SA));

median_error_margin_SA= [bootstrap_error_median_SA(ceil(alpha/2*B)),bootstrap_error_median_SA(ceil((1-alpha/2)*B))];
pos_frac_error_margin_SA= [bootstrap_error_pos_frac_gain_SA(ceil(alpha/2*B)),bootstrap_error_pos_frac_gain_SA(ceil((1-alpha/2)*B))];

display(['SA : gain on expected mse = ',num2str((1-(mean(mse_SA))/mean(mse_reg))*100,'%.1f%%'),' (+/- ',num2str((mse_error_margin_SA)/(1-mse_error_margin_reg)*100,'%.1f%%)') , ', average gain = ' num2str(mean(gains_SA)*100,'%.1f%%'),' (+/- ',num2str(mean_gain_error_margin_SA*100,'%.1f%%)'), ', median gain = ' num2str(mean((median(gains_SA)+median_error_margin_SA)*100),'%.1f%%'),' (+/- ',num2str(mean(abs(median_error_margin_SA))*100,'%.1f%%)'), ', fraction of positive gain = ' num2str(mean((mean(gains_SA>=0)+pos_frac_error_margin_SA)*100),'%.1f%%'),' (+/- ',num2str(mean(abs(pos_frac_error_margin_SA))*100,'%.1f%%)')]);
if visualization_1
figure;    
histogram(gains_SA,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 1.5]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with SA algorithm')
end
if visualization_2
figure;    
histogram((mean(mse_reg)-mse_SA)/mean(mse_reg),200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('gain vs expected MSE$_{reg}$, $g$','interpreter','latex')
%title('Histogram of the gain vs expected MSE$_{reg}$ obtained with SA algorithm','interpreter','latex')
end
end
if to_test(4)
mean_gain_error_margin_GA = norminv(1-alpha/2)*std(gains_GA)/(sqrt(n_draw_comp));
mse_error_margin_GA = norminv(1-alpha/2)*std(mse_GA)/(sqrt(n_draw_comp)*mean(mse_GA));

bootstrap_error_median_GA = sort(median(gains_GA )-bootstrp(B,@median,gains_GA ));
bootstrap_error_pos_frac_gain_GA = sort(mean(gains_GA >=0)-bootstrp(B,@(x) mean(x>=0),gains_GA ));

median_error_margin_GA = [bootstrap_error_median_GA(ceil(alpha/2*B)),bootstrap_error_median_GA(ceil((1-alpha/2)*B))];
pos_frac_error_margin_GA = [bootstrap_error_pos_frac_gain_GA(ceil(alpha/2*B)),bootstrap_error_pos_frac_gain_GA(ceil((1-alpha/2)*B))];

display(['GA : gain on expected mse = ',num2str((1-(mean(mse_GA ))/mean(mse_reg))*100,'%.1f%%'),' (+/- ',num2str((mse_error_margin_GA )/(1-mse_error_margin_reg)*100,'%.1f%%)') , ', average gain = ' num2str(mean(gains_GA )*100,'%.1f%%'),' (+/- ',num2str(mean_gain_error_margin_GA*100,'%.1f%%)'), ', median gain = ' num2str(mean((median(gains_GA )+median_error_margin_GA )*100),'%.1f%%'),' (+/- ',num2str(mean(abs(median_error_margin_GA ))*100,'%.1f%%)'), ', fraction of positive gain = ' num2str(mean((mean(gains_GA >=0)+pos_frac_error_margin_GA )*100),'%.1f%%'),' (+/- ',num2str(mean(abs(pos_frac_error_margin_GA ))*100,'%.1f%%)')]);

if visualization_1
figure;    
histogram(gains_GA,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 1.5]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA algorithm')
end
if visualization_2
figure;    
histogram((mean(mse_reg)-mse_GA)/mean(mse_reg),200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('gain vs expected MSE$_{reg}$, $g$','interpreter','latex')
%title('Histogram of the gain vs expected MSE$_{reg}$ obtained with GA algorithm','interpreter','latex')
end
end
if to_test(5)
mean_gain_error_margin_RT = norminv(1-alpha/2)*std(gains_RT)/(sqrt(n_draw_comp));
mse_error_margin_RT = norminv(1-alpha/2)*std(mse_RT)/(sqrt(n_draw_comp)*mean(mse_RT));
bootstrap_error_median_RT = sort(median(gains_RT )-bootstrp(B,@median,gains_RT ));
bootstrap_error_pos_frac_gain_RT = sort(mean(gains_RT >=0)-bootstrp(B,@(x) mean(x>=0),gains_RT ));

median_error_margin_RT = [bootstrap_error_median_RT(ceil(alpha/2*B)),bootstrap_error_median_RT(ceil((1-alpha/2)*B))];
pos_frac_error_margin_RT = [bootstrap_error_pos_frac_gain_RT(ceil(alpha/2*B)),bootstrap_error_pos_frac_gain_RT(ceil((1-alpha/2)*B))];

display(['RT : gain on expected mse = ',num2str((1-(mean(mse_RT ))/mean(mse_reg))*100,'%.1f%%'),' (+/- ',num2str((mse_error_margin_RT )/(1-mse_error_margin_reg)*100,'%.1f%%)') , ', average gain = ' num2str(mean(gains_RT )*100,'%.1f%%'),' (+/- ',num2str(mean_gain_error_margin_RT*100,'%.1f%%)'), ', median gain = ' num2str(mean((median(gains_RT )+median_error_margin_RT )*100),'%.1f%%'),' (+/- ',num2str(mean(abs(median_error_margin_RT ))*100,'%.1f%%)'), ', fraction of positive gain = ' num2str(mean((mean(gains_RT >=0)+pos_frac_error_margin_RT )*100),'%.1f%%'),' (+/- ',num2str(mean(abs(pos_frac_error_margin_RT ))*100,'%.1f%%)')]);

if visualization_1 
figure;
histogram(gains_RT,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 1.5]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with RT algorithm')
end
if visualization_2
figure;    
histogram((mean(mse_reg)-mse_RT)/mean(mse_reg),200,'Normalization','pdf')
xlim([-2, 1]);
ylim([0, 2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('gain vs expected MSE$_{reg}$, $g$','interpreter','latex')
%title('Histogram of the gain vs expected MSE$_{reg}$ obtained with RT algorithm','interpreter','latex')
end
end




%%
%mse_reg = 31.076762302039494; 
if sum(to_test/5)==1
fontsize = 7*2.5;
FIG = figure; hold on;
genList = (1:length(avgCostHist_GA));
evalList_RT = (1:n_eval);
evalList_GF = cumsum((GF_eval_max:-1:GF_eval_min)*n_eval_GF);
evalList_GF = round(evalList_GF*(n_eval/evalList_GF(end)));
evalList_GB = cumsum((GB_eval_max:-1:GB_eval_min)*n_eval_GB);
evalList_GB = round(evalList_GB*(n_eval/evalList_GB(end)));
the_end = round(genList*pop_size) ;
plot(round((genList-1)*pop_size),0*avgCostHist_GA+mean(mse_reg),':k');
plot(round(genList*pop_size),avgCostHist_GA,'--r');
plot(round(genList*pop_size),minCostHist_GA,'-*r');
%plot(round(genList*pop_size)*n_draw,avgCostHist_SA,'--m');
%plot((0:length(minCostHist_SA)-1)*the_end(end)/(length(minCostHist_SA)-1),minCostHist_SA,'-m');
plot(round(genList*pop_size),minCostHist_SA,'.-','color',[0.5,0.5,0.5]);

%plot(evalList_RT*n_draw,avgCostHist_RT,'--k');
plot(evalList_RT,minCostHist_RT(1:2:end),'-k');


plot(round(evalList_GF(end)*Greedy_plot_spacing),minCostEnd_GF(:),'^b');
plot(round(evalList_GB(end)*Greedy_plot_spacing),minCostEnd_GB(:),'v','color', [ 0.9100 0.4100 0.1700]);

box on
xlim([0 n_eval])
xticks([0:500:n_eval])
%ylim([1*floor(minCostHist_GA(end)/1) 1*ceil(minCostHist_RT(1)/1)])
%legend('regualar MSE cost','GA: average cost', 'GA: min cost', 'SA: min cost','RT: min cost','GF: min cost')
%legend({'regular cost', 'GA: min cost', 'SA: min cost','RT: min cost','GF: min cost','GB: min cost'},'interpreter','latex')
legend({'regualar cost','GA: average cost', 'GA: min cost', 'SA: min cost','RT: min cost','GF: min cost','GB: min cost'},'interpreter','latex')
%legend('regualar MSE cost','GA: average cost', 'GA: min cost', 'SA: min cost','RT: average cost','RT: min cost','GF: min cost','GB: min cost')
xlabel('number of cost function evaluation','interpreter','latex')
%ylabel('$\hat{\mathrm{E}}_{\mathrm{MSE}}$','interpreter','latex')
ylabel('Cost','interpreter','latex')

%title('Evolution of minimum cost','interpreter','latex')
%title('Evolution of the average and minimum cost $\hat{\mathrm{E}}_{\mathrm{MSE}}$','interpreter','latex')
set(findall(gcf,'-property','FontSize'),'Fontsize',fontsize);
set(0, 'DefaultLineLineWidth', 1);
end

%%
if to_test(4)
disp('Comparison over a single draw : GA and ref ');
meas_reg = round(linspace(0,T,n_measurements));

measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;

measurements_GA = zeros(1,T+1); 
measurements_GA(meas_GA+1) = 1;
figure
x0_single = initialization(1);
x_single = model(T,x0_single);
y_single = measurements(x_single,T,0);
part = initialization(n_part_fine,0,0);
tau_reg_single = particle_filter(y_single,measurements_reg,T,part,0);
tau_GA_single = particle_filter(y_single,measurements_GA,T,part,0);
err_reg_single = mean((objective(x_single,0)-tau_reg_single).^2,2);
err_GA_single = mean((objective(x_single,0)-tau_GA_single).^2,2);
 
gain_GA_single = (err_reg_single - err_GA_single)/err_reg_single;
display(['GA : gain = ' num2str((gain_GA_single)*100,'%.1f %%')]);

plot(0:0.25:T*0.25,(objective(x_single,0)),'k')
hold on
% plot(0:T,objective(tau_j_reg,0),':b')
% plot(0:T,objective(tau_j_GF,0),'-.r')
plot(0:0.25:T*0.25,tau_reg_single,':b')
plot(0:0.25:T*0.25,tau_GA_single,'-.r')
ylim([5*floor(min([min(tau_GA_single),min(objective(x_single,0)),min(tau_reg_single)])/5),5*ceil(max([max(tau_GA_single),max(objective(x_single,0)),max(tau_reg_single)])/5)])
y_limits_plot = ylim;
plot(meas_reg*0.25,y_limits_plot(2)*ones(n_measurements,1),'+b')
plot(meas_GA*0.25,y_limits_plot(1)*ones(n_measurements,1),'*r')
legend({'truth','$M_{reg}$','$M_{GA}$'},'interpreter','latex')
ylabel('z(t)','interpreter','latex')
xlabel('time, s','interpreter','latex')
xlim([0 T*0.25])
set(0,'defaultAxesFontSize',fontsize)
set(findall(gcf,'-property','FontSize'),'Fontsize',fontsize);


end
%%
if 1
T = 60;
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
plot(linspace(0,T/4,T+1),(objective(x_j,0)),'k')
ylabel('motion [mm] ','interpreter','latex')
xlim([0,T/4])
subplot(4,1,3);
hold on
ylabel('amplitude [mm] ','interpreter','latex')
plot(linspace(0,T/4,T+1),x_j(1,:),'k')
plot(linspace(0,T/4,T+1),a_minus*ones(1,T+1))
plot(linspace(0,T/4,T+1),a_plus*ones(1,T+1))
xlim([0,T/4])
subplot(4,1,4);
hold on
ylabel('offset [mm] ','interpreter','latex')
plot(linspace(0,T/4,T+1),b_minus*ones(1,T+1))
plot(linspace(0,T/4,T+1),b_plus*ones(1,T+1))
plot(linspace(0,T/4,T+1),x_j(3,:),'k')
xlim([0,T/4])


xlabel('time, [s]','interpreter','latex')
set(findall(gcf,'-property','FontSize'),'Fontsize',fontsize);
end