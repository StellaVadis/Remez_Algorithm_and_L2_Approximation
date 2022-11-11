warning('off')
%% Case 1
format long
clear
f = @(x) 2*x.^3+x.^2+x-1;
a = -1;
b = 1;
n = 4;
k = 0.0000001;
reference = [-1,-0.9,-0.8,-0.7,-0.6,-0.5];
threshold = 10^(-5);
epoch_max = 10000;
[p,p_table,reference_table,h_table] = Minimax_Approximation_Continuous(f,a,b,n,k,reference,threshold,epoch_max);
p

%% Case 2
format long
clear
f = @(x) 2*x.^3+x.^2+2*x-1;
a = -1;
b = 1;
n = 2;
k = 0.0000001;
reference = [-1,-0.9,-0.8,-0.7];
threshold = 10^(-5);
epoch_max = 10000;
[p,p_table,reference_table,h_table] = Minimax_Approximation_Continuous(f,a,b,n,k,reference,threshold,epoch_max);
p
for graphs = 1: size(reference_table,1)
    figure(graphs)
    plot(reference_table(graphs,:),f(reference_table(graphs,:)),'rs','Markersize',10)
    hold on
    fplot(f,[-1,1],'b')
    fplot(@(x) p_table(graphs,:) * x.^(0:size(p_table,2)-1)',[-1,1],'k')
    handles = legend('Reference Points','Approximation Function: f','Approximation Polynomial: p','location','northwest');
    set(handles,'FontName','Times New Roman','FontSize',10,'FontWeight','normal')
end

%% Case 3
format long
clear
f = @(x) exp(x);
a = -1;
b = 1;
n = 2;
k = 0.0000001;
reference = [-1,-0.9,-0.8,-0.7];
threshold = 10^(-11);
epoch_max = 10000;
[p,p_table,reference_table,h_table] = Minimax_Approximation_Continuous(f,a,b,n,k,reference,threshold,epoch_max);
p