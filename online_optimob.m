function [meas_elite,meas_1] = online_optimob(y,n_measurements,T,pop_size,max_gen,n_part,n_draw,measurements_spacing)

meas_elite = ones(1,n_measurements);
%For loop on the measurements to make, at each step we want to
%determine the following best measurement times, take a measure at the
%first one and iterate
% r = (T-floor(T/n_measurements))/T;
% % ratio = (1-r^(n_measurements))/(1-r);
% % r = 1; 
% ratio = 1; 
% n_part_online = round(n_part/(ratio)^(1/3));
% n_draw_online = round(n_draw/(ratio)^(1/3));
% pop_size_online = 2*round(pop_size*r/2);
% max_gen_online = round(max_gen/ratio^(1/3));

%[meas_1,~,~,~] = greedy_forward_algo(n_measurements,T,n_part,n_draw,measurements_spacing);
[meas_1,~,~,~] = genetical_algo(n_measurements,T,pop_size,max_gen,n_part,n_draw,measurements_spacing);
meas_elite(1) = meas_1(1);
meas_j_T = meas_1(2:end)-meas_1(1);


for j = 2:n_measurements
    T_j = T-meas_elite(j-1);
    n_measurements_j = n_measurements-(j-1);
    meas_1_j = meas_elite(1:(j-1));
    %[meas_j,~,~,~] = greedy_forward_algo(n_measurements_j,T_j,n_part,n_draw,measurements_spacing,y,meas_1_j);
    [meas_j,~,~,~] = genetical_algo(n_measurements_j,T_j,pop_size,max_gen,n_part,n_draw,measurements_spacing,y,meas_1_j,meas_j_T);
    meas_elite(j)=meas_elite(j-1)+meas_j(1);
    meas_j_T = meas_j(2:end)-meas_j(1);
%     pop_size_online = 2*round(pop_size_online*r/2);
end

end