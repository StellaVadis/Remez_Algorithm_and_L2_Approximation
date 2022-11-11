warning('off')
%% Case 1: 2-norm inner product
f = @(x) sqrt(1+x.^2);
w = @(x) 1;
a = 0;
b = 1;
n = 1;
threshold = 10^(-15);
integral_division = 100000;
[p,polynomial_table,phi_polynomial,alpha,beta] = Approximation_2_norm(f,w,a,b,n,threshold,integral_division);
p

%% Case 2: 2-norm Inner Product
f = @(x) exp(x);
w = @(x) 1;
a = -1;
b = 1;
n = 3;
threshold = 10^(-15);
integral_division = 100000;
[p,polynomial_table,phi_polynomial,alpha,beta] = Approximation_2_norm(f,w,a,b,n,threshold,integral_division);
p
for graphs = 1: n+1
    figure(graphs+10)
    fplot(f,[-1,1],'b')
    hold on
    fplot(@(x) polynomial_table(graphs,:) * x.^(0:size(polynomial_table,2)-1)',[-1,1],'k')
    handles = legend('Approximation Function: f','Approximation Polynomial: p','location','northwest');
    set(handles,'FontName','Times New Roman','FontSize',10,'FontWeight','normal')
end

%% Case 3: Weighted Inner Product, a more general case
f = @(x) exp(x);
w = @(x) abs(x);
a = -1;
b = 1;
n = 3;
threshold = 10^(-18);
integral_division = 100000;
[p,polynomial_table,phi_polynomial,alpha,beta] = Approximation_2_norm(f,w,a,b,n,threshold,integral_division);
p
for graphs = 1: n+1
    figure(graphs+20)
    fplot(f,[-1,1],'b')
    hold on
    fplot(@(x) polynomial_table(graphs,:) * x.^(0:size(polynomial_table,2)-1)',[-1,1],'k')
    handles = legend('Approximation Function: f','Approximation Polynomial: p','location','northwest');
    set(handles,'FontName','Times New Roman','FontSize',10,'FontWeight','normal')
end