%Newton Method
clear all;
syms x;

exp_f(x) = -exp(-x^2);
g_exp_f(x) = -2*(exp(-x^2))*x;
hess_exp_f(x) = exp(-x^2)*(2-4*(x^2));

% Point 1: x=0.1
k = 0;
k_max = 5;
status = 0;
x_k = 1;

%plot initial point
x_init = x_k;
y_init = exp_f(x_k);
hold on
fplot (exp_f(x),[-3 3])
plot(x_init,y_init,'b*','markersize',20)

while(status==0)
    old_x_k = x_k;
    d_k = - (g_exp_f(old_x_k))/(hess_exp_f(old_x_k));
    x_k = old_x_k + d_k;
    k = k + 1;
    
    fprintf ("\nx = %f",double (x_k));
    %plot each iteration point
    plot(x_k,exp_f(x_k),'ko-');
    
    %if(abs(hess_exp_f(x_k))<10e-3)
        %status =  -1;
    %end
    if(k==k_max)
        status =  -1;
    end
end

x_opt = double(x_k)
y_opt =  exp_f(x_k);
%plot the final point
plot(x_opt,y_opt,'r*','markersize',20)
hold off
