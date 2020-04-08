
T = 60;
n_measurements = 21;
n_draw_comp = 10000;
x0 = initialization(1,0,0);

n_draw = 100;
n_part = 200;

n_part_fine = 10000;

fontsize = 7*2;

%1.Simulation
%1.1. Random motion model : X_j
x = model(T,x0,0);
%1.2. Artificial data record Y_j
y = measurements(x,T,0);
meas_reg = round(linspace(0,T,n_measurements));


measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;

part = initialization(n_part_fine,0,0);
tau_reg = particle_filter(y,measurements_reg,T,part,0);
err_reg = mean((objective(x)-objective(tau_reg)).^2,2);

mse = zeros(1,n_draw_comp); 
display(['2. Computation MSE convergence by running ',num2str(n_draw_comp,'%.0f draws')]);
display('Computation computations started');
t_start_computation = tic;
for j = 1:n_draw_comp

    mse(j) = MC_MSE_estimator(meas_reg,T,n_draw,n_part,y,0);

end
t_elapsed_computation = toc(t_start_computation);
display(['Comparison computations completed, time elapsed = ',num2str(t_elapsed_computation,'%.0f sec')]);
%%
figure;
histogram((mse-mean(mse))/mean(mse),'Normalization','pdf')
%xlim([-2, 1]);
%ylim([0, 2.2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('normalized mse, $\hat{\mathrm{E}}_{\mathrm{MSE}}$','interpreter','latex')
title('Histogram of the normalized $\mathrm{MSE}$ estimator','interpreter','latex')
set(findall(gcf,'-property','FontSize'),'Fontsize',fontsize);

norm_mse = (mse-mean(mse))/mean(mse);