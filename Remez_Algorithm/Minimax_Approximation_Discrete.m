function [p,p_table,indices_table,reference_table,h_table,extremum_table] = Minimax_Approximation_Discrete(x,y,n,indices,threshold,epoch_max)

%% Input
% x: Abscissa of discrete points
% y: Ordinate of discrete points
% n: Order of Approximation polynomial
% indices: The initial indices of iterations
% threshold: When the inf-norm Approximation Error is less than threshold, the program terminates
% epoch_max: When the maximum epoch is achieved, the program terminates

%% Output
% p: The final approximation polynomial
% p_table: The approximation polynomial in each iteration
% indices_table: The indices of reference in each iteration
% reference_table: The reference in each iteration
% h_table: The levelled error h in each iteration
% extremum_table: The extremum of error function in each iteration

%% Handle the basic exception
if length(x) ~= length(y)
    error('The size of x and y is not the same!')
end

if length(x) < 0
    error('Invalid size x!')
end

if length(y) < 0
    error('Invalid size y!')
end
    
if n < 0
    error('Invalid n!')
end

if length(indices) ~= n+2
    error('Invalid size of initial indices and reference!')
end

if threshold < 0
    error('Invalid threshold value, please choose a non-zero threshold value!')
end

%%
for epoch = 1:epoch_max
    
    %% Establish the indices and reference Matrices
    reference = x(indices);
    indices_table(epoch,1:n+2) = indices;
    reference_table(epoch,1:n+2) = reference;
    
    %% Establish the coefficients Matrices A and beta
    A = ones(n+2,2);
    for i = 1:n+2
        A(i,1) = (-1)^i;
    end
    for j = 3:n+2
        for i = 1:n+2
            A(:,j) = reference.^(j-2);
        end
    end
    beta = y(indices)';
    A_beta = [A,beta];
    
    %% Gaussian Elimination to solve the simultaneous equations
    for centre = 1:n+2
        A_beta(centre,:) = A_beta(centre,:)./A_beta(centre,centre);
        for j = 1:size(A_beta,1)
            if centre~=j
                A_beta(j,:) = A_beta(j,:) - A_beta(centre,:)  *A_beta(j,centre);
            end
        end
    end
    
    %% The first end of row is the value of h, and the remaining are coefficeints of p in an ascending order(p_0,p_1,...,p_n)
    h = A_beta(1,end);
    h_table(epoch) = h;
    p = A_beta(2:end,end);
    p_table(epoch,1:n+1) = p';
    
    %% Find the x_enter with the largest magnitude of error under the p
    extremum_abs = -inf;
    key = 0;
    for i = 1:length(x)
        if abs(y(i)-x(i).^(0:n)* p) > extremum_abs
            extremum_abs = abs(y(i)-x(i).^(0:n)* p);
            extremum_sign = sign(y(i)-x(i).^(0:n)* p);
            key = i;
        end
    end
    extremum_table(epoch) = extremum_abs * extremum_sign;
    
    %% Termination Condition
    if (extremum_abs - abs(h)) < threshold
        disp('Terminates! We have reached the threshold!')
        break
    end
    
    %% Change the Reference
    if key < indices(1) && extremum_sign == -1*sign(h)
        indices(1) = key;
    elseif key < indices(1) && extremum_sign == sign(h)
        indices(2:n+2) = indices(1:n+1);
        indices(1) = key;
    elseif key > indices(n+2) && extremum_sign == (-1)^(n+2)*sign(h)
        indices(n+2) = key;
    elseif key > indices(n+2) && extremum_sign == (-1)^(n+1)*sign(h)
        indices(1:n+1) = indices(2:n+2);
        indices(n+2) = key;
    elseif key < indices(n+2) && key > indices(1)
        %To discover the location of x_enter
        %[1,2] means the x corresponding to key is located between references 1 and 2
        location = [1,2];
        while key > indices(location(2))
            location = location + 1;
        end
        
        if extremum_sign == (-1)^location(1)*sign(h)
            indices(location(1)) = key;
        elseif extremum_sign == (-1)^location(2)*sign(h)
            indices(location(2)) = key;
        end
    end
      
end