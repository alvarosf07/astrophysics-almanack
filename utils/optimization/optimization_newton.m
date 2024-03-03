%*******************************************************************************
%This Matlab function "Newton_method" implements the Newton method
%                                                                              *
%*****************************************************************************0*
%                                                    %**************************
%                                                    % Input Parameters        *
%                                                    %**************************
%
%            obj_f                % The objective function to minimize         %
%            starting_pt          % The starting point                         %
%            grad_f               % The gradient of the objective function     %
%            Hess_f               % The Hessian of the objective function      %
%            alpha_0              % The initial step-size                      %
%            k_max                % The maximum number of iterations           %
%            eps_tol              % The stopping convergence tolerance         %
%
%
%                                                    %**************************
%                                                    % Output Parameters       %
%                                                    %**************************
%
%            x_opt                % The optimal solution                       %
%            f_opt                % obj_f(x_opt)                               %
%            status               % An integer that indicates the termination
%                                   status of the method. I.e,
%                                   status =1  : the method converged
%                                   statut =-1 : the maximum number of iterations is reached
%            nit                  % The number of iterations performed by the  method
%            f_count              % The number of evaluation of obj_f          %
%            g_count              % The number of evaluation of grad_f         %
%
% Contact: Y. Diouane (youssef.diouane@isae.fr) --
% (C) Institut Sup�rieur de l'A�ronautique et de l'Espace (ISAE-Supa�ro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [                                                              ...
    x_opt  ,                                                            ...
    f_opt  ,                                                            ...
    status ,                                                            ...
    nit                                                                 ...
    ]                                                                   ...
    = Newton_method(                                                  ...
    obj_f           ,                                                   ...
    starting_pt     ,                                                   ...
    grad_f          ,                                                   ...
    Hess_f          ,                                                   ...
    k_max           ,                                                   ...
    eps_tol                                                             ...
    )

%**********************
%    CODE             *
%**********************

global f_count  ;                           % The number of evaluation of obj_f
global g_count  ;                           % The number of evaluation of grad_f
global h_count  ;                           % The number of evaluation of Hess_f

f_count    = 0                                                                 ;
g_count    = 0                                                                ;
h_count    = 0                                                                 ;

% Initialization
x_k        = starting_pt                                                       ;
n          = length( x_k  )                                                    ;
k          =   0                                                               ;
status     =   0                                                               ;
while(status==0)
    old_x_k  =  x_k                                                            ;
    %
    %  TO DO - COMPLETE THE UPDATE OF x_k
    x_k
    d_k      =  -inv(Hess_f(x_k))*grad_f(x_k)
    x_k      =  old_x_k + d_k
    k        =  k + 1
    
    %plot the current point
    plot([old_x_k(1) x_k(1) ], [old_x_k(2) x_k(2)],'ko-')                      ;
    
    %
    % TO DO - INCLUDE A STOPPING CRITERION
    %
    if (  norm(grad_f(x_k))<=eps_tol   )
        status =  -1                                                           ;
    end
    if(k==k_max)
        status=  -1                                                            ;
    end
end
x_opt    =   x_k                                                               ;
f_opt    =   feval(obj_f,x_opt)                                                ;
nit      =   k                                                                 ;
%plot the final point
plot(x_opt(1),x_opt(2),'b*','markersize',20)
return                                                                         ;
