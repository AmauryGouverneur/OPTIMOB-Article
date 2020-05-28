function part = mutation(part,t)
% The function "mutation" mutates the particles according to a discrete 
% stochastic nonlinear dynamical systems. 
% 
% x(t+1) = f(x(t),w(t)) for t = 0,...,T ? 1,
% 
% Inputs: 
%   - part : particles
%   - t : time 
% 
% Output: 
%   - part : mutated particles 
% 
% Implemented example: 
%   x(t+1) = x(t)/2 + 25*x(t)/(1+x(t)^2) + 8*cos(1.2t) + w(t)
%   w(t) ~ N(0,1)
%
% Date : 30/01/20
% Author : Amaury Gouverneur & Antoine Aspeel

n_part = size(part,2);
% part = part/2 + 25*part./(1+part.^2) + 8*cos(1.2*(t)) + randn(1,n_part)*sqrt(1);


a_minus = 8.8;
a_plus = 24;
sigma_a = 1;

omega_minus = 0.13*2*2;
omega_plus = 0.21*2*2;

sigma_omega = 0.005;

b_minus = -5.8;
b_plus = 5.8;
sigma_b = 1;


minus = [a_minus;omega_minus;b_minus];
plus = [a_plus;omega_plus;b_plus];
%sigmas = [sigma_a;sigma_omega/(t+1);sigma_b];
sigmas = [sigma_a;sigma_omega;sigma_b];

noise =  randn(3,n_part).*sigmas;
part = clip(part+noise,minus,plus);

end