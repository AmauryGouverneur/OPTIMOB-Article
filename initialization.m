function part = initialization(n_draw,y,meas_1_j)
% Function initializing the particles according to the model 
% X_0 ~ Xsi(x_0) dx_0
% 
% Input : 
%   - n_part : number of particles
% 
% Output : 
%   - part : particles
% 
% Date : 23/01/20
% Author : Amaury Gouverneur & Antoine Aspeel
if nargin == 1
    part = randn(1,n_draw)*5 ; 
else
    if meas_1_j == 0 
        part = randn(1,n_draw)*5 ;
    else 
        t_j = meas_1_j(end);
        index_t_j = t_j + 1;
        y_j = y(1:index_t_j);
        measurement_times = zeros(1,t_j+1); 
        measurement_times(meas_1_j+1) = 1;
        part = initialization(n_draw);
        [~,part] = particle_filter(y_j,measurement_times,t_j,part,0)  ; 
    end
end

end