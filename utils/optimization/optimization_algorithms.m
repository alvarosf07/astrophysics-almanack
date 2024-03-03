%**************************************************************************
% This is a program that compares different optimizatiom methods:
%
%    algo: 1, 2, or 3 to choose whcih method to run, i.e.,
%       algo=1   : The descent gradient method with a fixed step-size
%       algo=2   : The descent gradient method with the optimal step-size
%       algo=3   : The Newton method
%
%    prob: 1 or 2 to choose one of the test problem.
%       prob=1   : The function cone f1
%       prob=2   : The function rosenbrock f2
%
% Student: Álvaro Sánchez Fernández -- All Rights Reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************************
%     GLOBAL UPDATES      *
%**************************

global f_count  ;                      % The number of evaluation of obj_f
global g_count  ;                      % The number of evaluation of grad_f
global h_count  ;                      % The number of evaluation of Hess_f

fprintf ("\n\n\n")
fprintf ("****************************************************************************************************************\n");
fprintf ("                                             Beginning of the Program                                           \n");
fprintf ("****************************************************************************************************************\n");

%% Question 1: Computation of the Stationary Points
syms x y

fprintf ("\n")
fprintf ("               ************************************************************************************\n")
fprintf ("                                                   Question 1    \n");
fprintf ("               ************************************************************************************\n\n")

% Cone
eqn_cone = cone_2D (x,y); % We call the function to calculate the cone equation
gradient_cone = g_cone([x,y]); % We call the function to calculate the cone gradient
sol_cone = solve ([gradient_cone(1)==0,gradient_cone(2)==0],[x,y]); % Solve the system for x & y
stationary_points_cone = [sol_cone.x; sol_cone.y]; % Save the solution of the system in a vector
fprintf (" —> Stationary Points of The Cone Function: \n");
fprintf("     x =  %.2f \n", stationary_points_cone(1));
fprintf("     y =  %.2f \n", stationary_points_cone(2));

% Rosenbrock
eqn_rosenbrock = rosenbrock_2D (x,y); % We call the function to calculate the rosenberg equation
gradient_rosenbrock = g_rosenbrock([x,y]); % We call the function to calculate the rosenbrock gradient
sol_rosenbrock = solve ([gradient_rosenbrock(1)==0,gradient_rosenbrock(2)==0],[x,y]); % Solve the system for x & y
stationary_points_rosenbrock = [sol_rosenbrock.x,sol_rosenbrock.y]; % Save the solution of the system in a vector
fprintf ("\n —> Stationary Points of The Rosenbrock Function: \n");
fprintf("     x =  %.2f \n", stationary_points_rosenbrock(1));
fprintf("     y =  %.2f \n", stationary_points_rosenbrock(2));


% Question 2: Computation of the Stationary Points
fprintf ("\n\n")
fprintf ("               ************************************************************************************\n")
fprintf ("                                                   Question 2    \n");
fprintf ("               ************************************************************************************\n\n")

% Cone
hess_cone = h_cone ([stationary_points_cone(1),stationary_points_cone(2)]); % We call the function to calculate the cone hessian
eigenvalues_cone = eig (hess_cone);
fprintf ("    **************")
fprintf ("\n —> Cone Function: \n");
fprintf ("    **************\n")
fprintf ("\n     Eigenvalues of the Hessian Cone function : \n");
fprintf("     h1 =  %.0f \n", eigenvalues_cone(1));
fprintf("     h2 =  %.0f \n", eigenvalues_cone(2));
fprintf ("                    1) Condition 1 for Local Minimum: %.0f x %.0f > 0 \n",eigenvalues_cone(1),eigenvalues_cone(2));
fprintf ("                    2) Condition 2 for Local Minimum: %.0f + %.0f > 0 \n",eigenvalues_cone(1),eigenvalues_cone(2));

if det(hess_cone) > 0 && (trace(hess_cone) > 0)
    fprintf ("\n     Both Conditions are Satisfied: Point (x,y) = (%.0f,%.0f) is a Local Minimum for the Cone Function\n\n",stationary_points_cone(1),stationary_points_cone(2));
else
    fprintf ("\n     Point (x,y) = (%.0f,%.0f) is NOT a Local Minimum for the Cone Function\n\n",stationary_points_cone(1),stationary_points_cone(2));
end

% Rosenbrock
hess_rosenbrock = h_rosenbrock ([stationary_points_rosenbrock(1),stationary_points_rosenbrock(2)]); % We call the function to calculate the cone hessian
eigenvalues_rosenbrock = eig (hess_rosenbrock);
fprintf ("\n    ********************")
fprintf ("\n —> Rosenbrock Function: \n");
fprintf ("    ********************\n")
fprintf ("\n     Eigenvalues of the Hessian Rosenbrock function : \n");
fprintf("     h1 =  %.2f \n", eigenvalues_rosenbrock(1));
fprintf("     h2 =  %.2f \n", eigenvalues_rosenbrock(2));
fprintf ("                    1) Condition 1 for Local Minimum: %.2f x %.2f > 0 \n",eigenvalues_rosenbrock(1),eigenvalues_rosenbrock(2));
fprintf ("                    2) Condition 2 for Local Minimum: %.2f + %.2f > 0 \n",eigenvalues_rosenbrock(1),eigenvalues_rosenbrock(2));

if det(hess_rosenbrock) > 0 && (trace(hess_rosenbrock) > 0)
    fprintf ("\n     Both Conditions are Satisfied: Point (x,y) = (%.0f,%.0f) is a Local Minimum for the Rosenbrock Function\n\n",stationary_points_rosenbrock(1),stationary_points_rosenbrock(2));
else
    fprintf ("\n     Point (x,y) = (%.0f,%.0f) is NOT a Local Minimum for the Rosenbrock Function\n\n",stationary_points_rosenbrock(1),stationary_points_rosenbrock(2));
end

%% Question 5: Optimization Problem for Cone Function & Rosenbrock Function
fprintf ("\n")
fprintf ("\n")
fprintf ("               ************************************************************************************\n")
fprintf ("                                Question 5 - Gradient Descent Algorithm with Constant Step   \n");
fprintf ("               ************************************************************************************\n\n")


% Cone
fprintf ("\n    **************");
fprintf ("\n —> Cone Function: \n");
fprintf ("    **************\n");

algo = 1; % Algo 1 = Gradient Method, Fixed Step
prob = 1; % Problem 1 = Cone Function
Lab2_main (algo,prob);


fprintf ("\n\n\n    ********************");
fprintf ("\n —> Rosenbrock Function: \n");
fprintf ("    ******************** \n");

algo = 1; % Algo 1 = Gradient Method, Fixed Step
prob = 2; % Problem 2 = Rosenbrock Function
Lab2_main (algo,prob);


%% Question 8: Implement Gradient Descent Algorithm
fprintf ("\n")
fprintf ("\n")
fprintf ("               ************************************************************************************\n")
fprintf ("                               Question 8 - Gradient Descent Algorithm with Optimum Step   \n");
fprintf ("               ************************************************************************************\n\n")


% Cone
fprintf ("\n    **************");
fprintf ("\n —> Cone Function: \n");
fprintf ("    **************\n");

algo = 2; % Algo 2 = Gradient Method, Optimum Step
prob = 1; % Problem 1 = Cone Function
Lab2_main (algo,prob);


fprintf ("\n\n\n    ********************");
fprintf ("\n —> Rosenbrock Function: \n");
fprintf ("    ******************** \n");

algo = 2; % Algo 2 = Gradient Method, Optimum Step
prob = 2; % Problem 2 = Rosenbrock Function
Lab2_main (algo,prob);



%% Question 9: Implement the Newton Method
fprintf ("\n")
fprintf ("\n")
fprintf ("               ************************************************************************************\n")
fprintf ("                                         Question 9 - Newton Algorithm      \n");
fprintf ("               ************************************************************************************\n\n")


% Cone
fprintf ("\n    **************");
fprintf ("\n —> Cone Function: \n");
fprintf ("    **************\n");

algo = 3; % Algo 3 = Newton Method
prob = 1; % Problem 1 = Cone Function
Lab2_main (algo,prob);


fprintf ("\n\n\n    ********************");
fprintf ("\n —> Rosenbrock Function: \n");
fprintf ("    ******************** \n");

algo = 3; % Algo 3 = Newton Method
prob = 2; % Problem 2 = Rosenbrock Function
Lab2_main (algo,prob);


%% Question 9: Implement the Newton Method
fprintf ("\n")
fprintf ("\n")
fprintf ("               ************************************************************************************\n")
fprintf ("                                 Question 11 - Newton Algorithm for Exponential Function      \n");
fprintf ("               ************************************************************************************\n\n")


% Cone
fprintf ("\n  *********************");
fprintf ("\n —> Exponential Function: \n");
fprintf ("   *********************\n");

algo = 3; % Algo 3 = Newton Method
prob = 3; % Problem 3 = Exponential Function
Lab2_main (algo,prob);
