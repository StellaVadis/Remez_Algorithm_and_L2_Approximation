function [p,polynomial_table,phi_polynomial,alpha,beta] = Approximation_2_norm(f,weight,a,b,n,threshold,integral_division)

%% Input
% f: Approximation function
% weight : The weight function of inner product, if you want to use 2-norm approximation, use weight = @(x) 1
% a: Left boundary point of Approximation Interval [a,b]
% b: Right boundary point of Approximation Interval [a,b]
% n: Highest Order of Approximation Polynomial
% threshold: When the 2-norm Approximation Error is less than threshold, the program terminates
% integral_division: In Numerical Integration, the number of division of sub-intervals

%% Output
% p: The final approximation polynomial
% polynomial_table: The approximation polynomial in each epoch
% phi_polynomial: The phi_polynomial in each epoch
% alpha: The alpha in each epoch
% beta: The beta in each epoch

%% Handle the basic exception
if n < 0
    error('Invalid n!')
end

if threshold < 0
    error('Invalid threshold value, please choose a non-zero threshold value!')
end

if a > b
    error('Invalid interval, please choose a < b!')
end

if integral_division < 10000
    warning('The integral division is not dense enough, please choose an integral_division larger than 10000!')
end

%% If n == 0, the polynomial is order 0. Remark: phi{1} in the Matlab records phi_0, as we do not have 0-index in Matlab
%If n == 0, only this part of code will be executed
if n >= 0
    %This is used to control the termination condition
    break_console = 0;
    phi{1} = @(x) 1;
    phi_polynomial(1) = 1;
    phi_norm_square(1) = Numerical_Inner_Product(phi{1},phi{1},weight,a,b,integral_division);
    x_multiply_phi{1} = @(x) phi{1}(x)*x;
    alpha(1) = Numerical_Inner_Product(x_multiply_phi{1},phi{1},weight,a,b,integral_division)/phi_norm_square(1);
    c(1) = Numerical_Inner_Product(phi{1},weight,f,a,b,integral_division)/phi_norm_square(1);
    p(1,1:1) = c(1)*phi_polynomial(1,:);
    polynomial_table(1,1:1) = p;
    beta(1) = NaN;
    
    if abs(Numerical_Inner_Product(@(x) f(x) - c(1) * phi{1}(x),@(x) f(x) - c(1) * phi{1}(x),weight,a,b,integral_division)) < threshold
        Numerical_Inner_Product(@(x) f(x) - c(1) * phi{1}(x),@(x) f(x) - c(1) * phi{1}(x),weight,a,b,integral_division)
        break_console = 1;
        disp('Terminates! We have reached the threshold!')
    end
    
end

%% If n == 1, the polynomial is order 1
%If n == 1, only the codes before and this part will be executed
if n >= 1 && break_console ~= 1
    phi_polynomial(2,1:2) = 0;
    phi_polynomial(2,1:1) = phi_polynomial(2,1:1) - alpha(1) * phi_polynomial(1,1:1);
    phi_polynomial(2,2:2) = phi_polynomial(2,2:2) + phi_polynomial(1,1:1);
    phi{2} = @(x) phi_polynomial(2,1:2) * x.^(0:1)';
    phi_norm_square(2)= Numerical_Inner_Product(phi{2},phi{2},weight,a,b,integral_division);
    x_multiply_phi{2} = @(x) phi{2}(x)*x;
    alpha(2) = Numerical_Inner_Product(x_multiply_phi{2},phi{2},weight,a,b,integral_division)/phi_norm_square(2);
    beta(2) = phi_norm_square(2)/phi_norm_square(1);
    c(2) = Numerical_Inner_Product(phi{2},f,weight,a,b,integral_division)/phi_norm_square(2);
    p(1,1:2) = [p,0] + c(2) * phi_polynomial(2,:);
    polynomial_table(2,1:2) = p;
    
    if abs(Numerical_Inner_Product(@(x) f(x) - p(1:2)*x.^(0:1)',@(x) f(x) - p(1:2)*x.^(0:1)',weight,a,b,integral_division)) < threshold
        break_console = 1;
        disp('Terminates! We have reached the threshold!')
    end
    
end


%% If n >= 2, all the codes will be executed
if n >= 2 && break_console ~= 1
    for j = 3:n + 1
        phi_polynomial(j,1:j) = 0;
        phi_polynomial(j,2:j) = phi_polynomial(j-1,1:j-1);
        phi_polynomial(j,1:j) = phi_polynomial(j,1:j) - phi_polynomial(j-1,1:j) * alpha(j-1);
        phi_polynomial(j,1:j) = phi_polynomial(j,1:j) -beta(j-1) * phi_polynomial(j-2,:);
        phi{j} = @(x) phi_polynomial(j,1:j) * x.^(0:j-1)';
        
        phi_norm_square(j)= Numerical_Inner_Product(phi{j},phi{j},weight,a,b,integral_division);
        x_multiply_phi{j} = @(x) phi{j}(x)*x;
        alpha(j) = Numerical_Inner_Product(x_multiply_phi{j},phi{j},weight,a,b,integral_division)/phi_norm_square(j);
        beta(j) = phi_norm_square(j)/phi_norm_square(j-1);
        c(j) = Numerical_Inner_Product(phi{j},f,weight,a,b,integral_division)/phi_norm_square(j);
        p(1,1:j) = [p,0] + c(j)*phi_polynomial(j,:);
           
        polynomial_table(j,1:j) = p;
        
        if abs(Numerical_Inner_Product(@(x) f(x) - p(1:j)*x.^(0:j-1)',@(x) f(x) - p(1:j)*x.^(0:j-1)',weight,a,b,integral_division)) < threshold
            disp('Terminates! We have reached the threshold!')
            break
        end
        
    end
end