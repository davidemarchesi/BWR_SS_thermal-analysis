function I = trap_int_powerpurp(a, b, N, fun)
% Integral Calculated with the Trapezium method
% Inputs:
% - a,b: Integration Extremes;
% - N: Number of subintervals;
% - fun: function;
% Output:
% - I: integral;

% Subintervals length
h = (b-a)/N;

% Subintervals extremes/medium points
x = a : h : b;

% Evaluation of the points
fx = fun(x);

% Approximated integral calculation
I = 0;
for j = 1:N
    I = I + h.* (fx(j)+fx(j+1)).*0.5;
end
end