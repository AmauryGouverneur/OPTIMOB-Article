function [meas_elite] = online_optimob(y,n_measurements,T,pop_size,max_gen,n_part,n_draw)

meas_elite = ones(1,n_measurements);
%For loop on the measurements to make, at each step we want to
%determine the following best measurement times, take a measure at the
%first one and iterate

[meas_1,~,~,~] = genetical_algo(n_measurements,T,pop_size,max_gen,n_part,n_draw);
meas_elite(1) = meas_1(1);

for j = 2:n_measurements
    T_j = T-meas_elite(j-1);
    n_measurements_j = n_measurements-(j-1);
    meas_1_j = meas_elite(1:(j-1));
    [meas_j,~,~,~] = genetical_algo(n_measurements_j,T_j,pop_size,max_gen,n_part,n_draw,y,meas_1_j);
    meas_elite(j)=meas_elite(j-1)+meas_j(1);
end

end