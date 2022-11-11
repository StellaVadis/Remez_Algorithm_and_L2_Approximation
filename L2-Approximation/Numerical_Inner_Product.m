function values = Numerical_Inner_Product(f1,f2,w,a,b,integral_division)

%% Input
% f1: The first function of inner product <f1,f2>
% f2: The second function of inner product <f1,f2>
% w: weight function 
% a: Left boundary point of integral interval [a,b]
% b: Right boundary point of integral interval [a,b]
% integral_division: In Numerical Integration, the number of division of sub-intervals

%% Output
% values: The weighted inner product

%% Integral and Inner Product by Riemann Integral
values = 0;
g = @(x) w(x)*f1(x)*f2(x);
for i = 1:integral_division
    values = values + (b-a)/integral_division/2*(g(a+(b-a)*i/integral_division)+g(a+(b-a)/integral_division*(i-1)));
end