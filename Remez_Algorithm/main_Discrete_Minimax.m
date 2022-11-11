warning('off')
%% Case 1
format short
clear
x = -1:0.01:1;
y = 2*x.^3+x.^2+x-1;
n = 4;
indices = 1:n+2;
threshold = 10^(-12);
epoch_max = 10000;
[p,p_table,indices_table,reference_table,h_table,extremum_table] = Minimax_Approximation_Discrete(x,y,n,indices,threshold,epoch_max);
p

%% Case 2
format short
clear
x = -1:0.1:1;
y = 2*x.^3+x.^2+2*x-1;
n = 2;
indices = 1:n+2;
threshold = 10^(-12);
epoch_max = 10000;
[p,p_table,indices_table,reference_table,h_table,extremum_table] = Minimax_Approximation_Discrete(x,y,n,indices,threshold,epoch_max);
p
for graphs = 1: size(reference_table,1)
    figure(graphs)
    not_indices = 1:length(x);
    not_indices(indices_table(graphs,:)) = [];
    plot(x(indices_table(graphs,:)),y(indices_table(graphs,:)),'rs','Markersize',10)
    hold on
    plot(x(not_indices),y(not_indices),'bs','Markersize',10)
    fplot(@(x) p_table(graphs,:) * x.^(0:n)',[-1,1],'k')
    handles = legend('Reference Points','Points except Reference Points','Approximation Polynomial: p','location','northwest');
    set(handles,'FontName','Times New Roman','FontSize',10,'FontWeight','normal')
end

%% Case 3
format long
clear
x = -1:0.1:1;
y = exp(x);
n = 2;
indices = 1:n+2;
threshold = 10^(-12);
epoch_max = 10000;
[p,p_table,indices_table,reference_table,h_table,extremum_table] = Minimax_Approximation_Discrete(x,y,n,indices,threshold,epoch_max);
p