function [ xvect, it ] = newton_method_mod( x0 , nmax , tol , f , df, m )
% 
% [ xvect, it ] = newton_method (x0 , nmax , tol , f , df )
%
% This function has the aim to use the newton's method to solve non linear
% equations f(x)=0 ;
%
% Input parameters:
% x0 : starting point for the analysis;
% nmax : maximum number of iterations admitted;
% tol : the tolerance we want for the solution;
% f : the function;
% df : the derivate of the function;
%
% Output parameters:
% xvect : iterations results (the last is the solution);
% it : number of iterations;

% Output and variables inizialization
it=0;
xvect=[];
deltax=1000;
res=tol+1;
xold = x0 + 1000;

while ( it<nmax && deltax>=tol && res>=tol )
    xold = x0;
    if abs(df(x0))< eps
        disp('Stop of the method due to df(xk)=0 \n');
        it = it + 1;
        break
    end
    x0 = x0 - ( m.* (f(x0)/df(x0)));
    xvect = [xvect ; x0];
    res = abs(f(x0));
    deltax = xold - x0;
    it = it + 1;
end

if (it<nmax)
    fprintf('Convergence at the iteration k : %d \n ', it);
else
    fprintf('Maximum number of iterations reached : it = %d = nmax \n ', nmax);
end

fprintf(' Value calculated : %-12.8f \n', xvect(it));

return

% This function can also be transformed into the modified-newton-method adding to the inputs a variable 'm'
% which represents the multiplicity of the function in the zero (that must be known in
% order to use it), and changing at the line 33 the equation present with:
% x0 = x0 - ( m.* (f(x0)/df(x0)));