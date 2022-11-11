function [p,p_table,reference_table,h_table] = Minimax_Approximation_Continuous(f,a,b,n,k,reference,threshold,epoch_max)

%% Input
% f: Approximation function
% a: Left boundary point of Approximation Interval [a,b]
% b: Right boundary point of Approximation Interval [a,b]
% n: Order of Approximation polynomial
% k: The mini step value in calculating the derivative e(x)
% reference: The initial reference of iteration
% threshold: When the inf-norm Approximation Error is less than threshold, the program terminates
% epoch_max: When the maximum epoch is achieved, the program terminates

%% Output
% p: The final approximation polynomial
% p_table: The approximation polynomial in each iteration
% reference_table: The reference in each iteration
% h_table: The levelled error h in each iteration

%% Handle the basic exception
if a > b
    error('Invalid interval, please choose a < b!')
end

if n < 0
    error('Invalid n!')
end

if k < 0
    error('Invalid k!')
end

if length(reference) ~= n+2
    error('Invalid size of initial reference!')
end

if threshold < 0
    error('Invalid threshold value, please choose a non-zero threshold value!')
end

%%
for epoch = 1:epoch_max
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
    A_beta = [A,f(reshape(reference,length(reference),1))];
    
    %% Gaussian Elimination to solve the simultaneous equations
    for centre = 1:n+2
        A_beta(centre,:) = A_beta(centre,:)./A_beta(centre,centre);
        for j = 1:size(A_beta,1)
            if centre~=j
                A_beta(j,:) = A_beta(j,:) - A_beta(centre,:)*A_beta(j,centre);
            end
        end
    end
    h = A_beta(1,end);
    p = A_beta(2:end,end);
    
    reference_table(epoch,1:n+2) = reference;
    h_table(epoch) = h;
    p_table(epoch, 1:n+1) = p;
    
    %% Termination Condition:  When the increase of |h| in each epoch is smaller than the threshold, the program terminates
    if epoch >= 2
        if h_table(epoch) - h_table(epoch -1) < threshold
            break
        end
    end
    
    %% Define the error function and derivative of error function
    e = @(x) f(x) - p'*x.^(0:n)';
    e_derivative = @(x) (e(x+k)-e(x-k))/2/k;
    
    %% The extremum sometimes occurs on the boundary points
    % The small value 0.0000001 is use to control the rounding error
    if abs(e(a)) > abs(h) + 0.0000001 && abs(e(a)) >= abs(e(b)) && sign(e(a)) == sign(e(reference(1)))
        reference(1) = a;
    elseif abs(e(a)) > abs(h)  + 0.0000001  && abs(e(a)) >= abs(e(b)) && sign(e(a)) == -sign(e(reference(1)))
        reference(2:end) = reference(1:end-1);
        reference(1) = a;
    elseif abs(e(b)) > abs(h)  + 0.0000001  && abs(e(b)) > abs(e(a)) && sign(e(b)) == sign(e(reference(end)))
        reference(end) = b;
    elseif abs(e(b)) > abs(h)  + 0.0000001  && abs(e(b)) > abs(e(a)) && sign(e(b)) == -sign(e(reference(end)))
        reference(1:end-1) = reference(2:end);
        reference(end) = b;
        
        %% Search for a local maximum or minimum of error function that larger than h or less than -h by continuous binary search
    else
        break_console = 0;
        for sub_epoch = 0:10
            basis = 0:1/2^sub_epoch:1;
            s = a + basis*(b-a);
            for i = 1:length(s)-1
                if sign(e_derivative(s(i))) == -sign(e_derivative(s(i+1))) && sign(e_derivative(s(i)))~=0 && sign(e_derivative(s(i+1)))~=0
                    left = s(i);
                    right = s(i+1);
                    for searching = 1:1000
                        if e_derivative((left+right)/2) == 0
                            %disp('well done')
                            x0 = (left+right)/2;
                            break
                        elseif sign(e_derivative((left+right)/2)) == sign(e_derivative(right))
                            right = (left+right)/2;
                        elseif sign(e_derivative((left+right)/2)) == sign(e_derivative(left))
                            left = (left+right)/2;
                        else
                            %disp('dwduhwodq')
                        end
                    end
                    x0 = (left+right)/2;
                    if abs(e(x0)) > abs(h)  + 0.0000001 %this small number is used to control the deviation of derivative
                        break_console = 1;
                        break
                    end
                end
            end
        end
        
        %% Change the reference
        if x0 < reference(1) && sign(e(x0)) == sign(e(reference(1)))
            reference(1) = x0
        elseif x0 < reference(1) && sign(e(x0)) == -sign(e(reference(1)))
            reference(2:end) = reference(1:end-1);
            reference(1) = x0
        elseif x0 > reference(end) && sign(e(x0)) == sign(e(reference(end)))
            reference(end) = x0
        elseif x0 > reference(end) && sign(e(x0)) == -sign(e(reference(1)))
            reference(1:end-1) = reference(2:end);
            reference(end) = x0
        elseif x0 > reference(1) && x0 < reference(end)
            %To discover the location of key
            %[1,2] means the x corresponding to key is located between references 1 and 2
            location = [1,2];
            while x0 > reference(location(2))
                location = location + 1;
            end
            
            if sign(e(x0)) == (-1)^location(1)*sign(h)
                reference(location(1)) = x0;
            elseif sign(e(x0)) == (-1)^location(2)*sign(h)
                reference(location(2)) = x0;
            end
            %The other cases may be one of e(x0) or sign(h) is zero, which
            %means we have achieved a zero-error approximation
        else
            disp('The excellent zero error is achieved!')
        end
    end
    
end




