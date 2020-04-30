function x = model_2(T,x0)
% The function "model" simulates a realisation of a state vector according 
% to a discrete stochastic nonlinear dynamical systems. 
% 
% x(0) ~ F
% x(t+1) = f(x(t),w(t)) for t = 0,...,T ? 1,
%
% where w_t is the process noise. 
%
% Input: 
%   - T : length of the time interval  
% 
% Output: 
%   - x : state vector 
%
% Implemented example:
%   x(t+1) = x(t)/2 + 25*x(t)/(1+x(t)^2) + 8*cos(1.2t) + w(t)
%   x(0) ~ N(0,5^2)
%   w(t) ~ N(0,1)
% 
% Date : 30/01/20
% Author : Amaury Gouverneur & Antoine Aspeel

a_minus = 8.8;
a_plus = 24;
sigma_a = 0.1;

omega_minus = 0.13;
omega_plus = 0.21;
sigma_omega = 0.01;

b_minus = -5.8;
b_plus = 5.8;
sigma_b = 0.1;

minus = [a_minus;omega_minus;b_minus];
plus = [a_plus;omega_plus;b_plus];
sigmas = [sigma_a;sigma_omega;sigma_b];

if nargin ==1 
    x0 = (minus+plus)/2 + randn(3,1).*sigmas;
end

x = zeros(3,T+1);
x(:,0+1) = x0;

for t = 0:T-1
    index_t = t+1;
    noise =  randn(3,1).*sigmas;
    x(:,index_t+1) = clip(x(:,index_t)+noise,minus,plus);
end

end

function res = clip(x,a,b)

res = x.*(x>a & x<b) + a.*(x<=a) + b.*(x>=b);
end