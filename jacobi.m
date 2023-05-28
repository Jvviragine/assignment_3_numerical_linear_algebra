function [x, n_iterations, xs] = jacobi(A, b, x1, E, order)
%Approximates the solution for Ax = b using the Jacobi Method up to a certain Tolerance

%INPUT:
% 'A' is the Coeficient Matrix
% 'b' is the right hand side Vector(Column Vector)
% 'x1'is the initial guess for the solution(Column Vector)
% 'E' is the tolerance
% 'order' is the order in which the variables for the solution will be
% updated. It can be '1' for in order or '-1' for revserved order. If
% another number is specified, reverse order will be performed (assignment)


% OUTPUT: 
% 'x' is the final approximation for the solution(Column Vector)
% 'n_iterations' is the number of approximations needed to get to the tolerance 'E' 
% 'xs' are all the approximations for 'x' generated, including 'x' itself

% 1- Checks if A, b, and x1 have the same Dimension (Columns)
if (size(A, 1) == length(b)) && (length(b) == length(x1)) && (~isempty(x1))

    xs = zeros(length(b), 100); % Allocating a certain Memory for the Previous X Approximations
    % I chose 100 because for most approximations, it will be more than
    % enough

    n_iterations = 0; % We know at this point we will do at least one iteration, so we allocate the memory

    % Defines variables to store the previous iteration and next
    % approximations
    x_previous = x1;
    x_next = zeros(length(b), 1);

    % Defines a Vector that will store the difference between x_next and
    % x_previous
    x_difference = zeros(length(b), 1);

    % Adds the first guess to the list 
    xs(:, 1) = x1;

    % Define the Stop conditions for the Iteration
    while (n_iterations == 0 || norm((x_difference), inf) > E)
        % While we have not done any iteration or we are above the Tolerance

        % Checks which order the user desires (In-order or Inverse)
        if order == 1 % In Order
            
            for i = 1:length(b) % Going through the 'b' values
                
                x_next(i) = b(i);
                % Go through the values from Aij*xj
                for j = 1:length(b)
             
                    % Checks if i != j (we do not subtract when =)
                    if i ~= j
                        x_next(i) = x_next(i) - (A(i, j)*x_previous(j));
                    end
                end
                % Divide x_next by aii
                x_next(i) = x_next(i) ./ A(i, i);
            end
            n_iterations = n_iterations + 1; % Iterations Counter

            x_difference = x_next - x_previous;

            x_previous = x_next;

            xs(:, n_iterations + 1) = x_next;

            x_next = zeros(length(b), 1); 

        else % This means that it's a Reverse Order

            for i = length(b):-1:1 % Going through the 'b' values
                
                x_next(i) = b(i);
                % Go through the values from Aij*xj
                for j = length(b):-1:1
             
                    % Checks if i != j (we do not subtract when =)
                    if i ~= j
                        x_next(i) = x_next(i) - (A(i, j)*x_previous(j));
                    end
                end
                % Divide x_next by aii
                x_next(i) = x_next(i) ./ A(i, i);
            end
            n_iterations = n_iterations + 1; % Iterations Counter

            x_difference = x_next - x_previous;

            x_previous = x_next;

            xs(:, n_iterations + 1) = x_next;

            x_next = zeros(length(b), 1); 

        end 

    end
    % At this point, the value has converged, so we can assign the final
    % values for the variables
    x = x_previous;

    % See how many Iterations we used and remove the unused values from
    % 'xs'
    xs = xs(:, 1:(n_iterations + 1));

        

else % The case where the dimensions do not match
    x = 0; n_iterations = 0; xs = 0;
    disp('Double check the dimensions of the Inputs!')

end 


end

