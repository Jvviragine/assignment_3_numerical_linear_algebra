% DRIVER SCRIPT FOR ASSIGNMENT 3: NUMERICAL LINEAR ALGEBRA
clear variables; clear figures; clc;
format longG;

%% Setting up the Variables for the Problem
D_values = 1:5;
% All the Variations for the A Matrix
A1 = [-11, D_values(1), -1; -2, 5, D_values(1); D_values(1), 3, 7];
A2 = [-11, D_values(2), -1; -2, 5, D_values(2); D_values(2), 3, 7];
A3 = [-11, D_values(3), -1; -2, 5, D_values(3); D_values(3), 3, 7];
A4 = [-11, D_values(4), -1; -2, 5, D_values(4); D_values(4), 3, 7];
A5 = [-11, D_values(5), -1; -2, 5, D_values(5); D_values(5), 3, 7];
% Setting up the b Vector(Column Vector)
b = [3; 4; 1];
% Setting up the Initial Guess (0 Vector)
x1 = [0; 0; 0];
% Setting up the Order (-1 for Inverse)
order = -1;
% Setting up the Tolerance
E = 10^-6;

%% Testing the Jacobi Method for all the Matrices using Inverse Order
% Storing the Approximations
[x_jacobi_A1, n_iterations_jacobi_A1, xs_jacobi_A1] = jacobi(A1, b, x1, E, -1);
[x_jacobi_A2, n_iterations_jacobi_A2, xs_jacobi_A2] = jacobi(A2, b, x1, E, -1);
[x_jacobi_A3, n_iterations_jacobi_A3, xs_jacobi_A3] = jacobi(A3, b, x1, E, -1);
[x_jacobi_A4, n_iterations_jacobi_A4, xs_jacobi_A4] = jacobi(A4, b, x1, E, -1);
[x_jacobi_A5, n_iterations_jacobi_A5, xs_jacobi_A5] = jacobi(A5, b, x1, E, -1);



%% Testing the Gauss-Seidel for all the Matrices using Inverse Order
% Storing the Approximations
[x_gauss_seidel_A1, n_iterations_gauss_seidel_A1, xs_gauss_seidel_A1] = gauss_seidel(A1, b, x1, E, -1);
[x_gauss_seidel_A2, n_iterations_gauss_seidel_A2, xs_gauss_seidel_A2] = gauss_seidel(A2, b, x1, E, -1);
[x_gauss_seidel_A3, n_iterations_gauss_seidel_A3, xs_gauss_seidel_A3] = gauss_seidel(A3, b, x1, E, -1);
[x_gauss_seidel_A4, n_iterations_gauss_seidel_A4, xs_gauss_seidel_A4] = gauss_seidel(A4, b, x1, E, -1);
[x_gauss_seidel_A5, n_iterations_gauss_seidel_A5, xs_gauss_seidel_A5] = gauss_seidel(A5, b, x1, E, -1);



%% Data for Display
x_jacobi_A1
n_iterations_jacobi_A1
x_gauss_seidel_A1
n_iterations_gauss_seidel_A1
x_exact_A1 = A1\b


x_jacobi_A2
n_iterations_jacobi_A2
x_gauss_seidel_A2
n_iterations_gauss_seidel_A2
x_exact_A2 = A2\b

x_jacobi_A3
n_iterations_jacobi_A3
x_gauss_seidel_A3
n_iterations_gauss_seidel_A3
x_exact_A3 = A3\b

x_jacobi_A4
n_iterations_jacobi_A4
x_gauss_seidel_A4
n_iterations_gauss_seidel_A4
x_exact_A4 = A4\b

x_jacobi_A5
n_iterations_jacobi_A5
x_gauss_seidel_A5
n_iterations_gauss_seidel_A5
x_exact_A5 = A5\b

%% Question 4 - Analysing and Plotting the Norm for d = 3

% Generating the Differences between subsequent Iterations for Jacobi
differences_jacobi_norm = zeros(1, n_iterations_jacobi_A3); % Allocating Space to store the Norms for Jacobi

% Calculating the Norm between Subsequent Approximations for Jacobi
for i = 1:(n_iterations_jacobi_A3)
    
    differences_jacobi_norm(i) = norm((xs_jacobi_A3(:, i) - xs_jacobi_A3(:, i+1)));
end

% Generating the Differences between subsequent Iterations for Gauss-Seidel
differences_gauss_seidel_norm = zeros(1, n_iterations_gauss_seidel_A3); % Allocating Space to store the Norms for Gauss-Seidel

% Calculating the Norm between Subsequent Approximations for Gauss-Seidel
for i = 1:(n_iterations_gauss_seidel_A3)
    
    differences_gauss_seidel_norm(i) = norm((xs_gauss_seidel_A3(:, i) - xs_gauss_seidel_A3(:, i+1)));
end

%% Plotting the Norms of the Differences for both Methods with d = 3

% Plotting for Jacobi
figure(1)
plot(1:n_iterations_jacobi_A3, differences_jacobi_norm, "--bo")
grid on
xlabel("Iteration"); ylabel("Two Norm of the Difference"); title("Norm of the Difference of Subsequent Iterations for Jacobi with d = 3");

% Plotting for Gauss-seidel
figure(2)
plot(1:n_iterations_gauss_seidel_A3, differences_gauss_seidel_norm, "--ko")
grid on
xlabel("Iteration"); ylabel("Two Norm of the Difference"); title("Norm of the Difference of Subsequent Iterations for Gauss-Seidel with d = 3");
