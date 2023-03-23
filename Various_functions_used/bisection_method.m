function [xvect,xdif,fx,it]=bisection_method(a,b,nmax,toll,f)
%
% [xvect,xdif,fx,it]=bisection_method(a,b,nmax,toll,f) 
%
% Bisection method for the resolution of the non linear equation f(x)=0
%
% Input parameters:
%
% a,b       Research range extremes
% nmax      Maximum number of requested iterations
% toll      Tolerance
% f         Function handle containing the function
%
% Output parameters:
%
% xvect     Vect. containing all the results from every iteration
%           (the last is the solution)
% xdif      Vect. containing the progress between two consecutive iterations
% fx        Vect. containing 'f' values of 'xvect'
% it        Iterations made


% output parameters inizialization
it=0; % the first step will be 0
xvect=[];
fx=[]; 
xdif=[];
delta_x=0;

% until the iteration is not started, i don't know the error, so in order
% to start the iteration a fictitious value for the error is set
err=toll+1; 


% First we verify if the zeros theorem is applicable into the interval
% extremes considerated for the resolution
if f(a)*f(b)>0
    disp("The function has no zero in the defined interval \n");
else
% The iteration continues until i don't reach the maximum number of
% iterations requested or i reach an error value under the tolerance
        while (it < nmax && err > toll)
            
            if it==1
                if f(a)==0
                    xvect=[xvect,a];
                    disp("The 'A' interval extreme is a zero");
                    break
                end
                if f(b)==0
                    xvect=[xvect,b];
                    disp("The 'B' interval extreme is a zero");
                    break
                end
            end
            
            x_m=(a+b)/2;
            
            if ( f(a)*f(x_m)<0 )
                delta_x=abs(b-x_m);
                b=x_m;
            end
            if ( f(b)*f(x_m)<0 )
                delta_x=abs(a-x_m);
                a=x_m;
            end
            xvect=[xvect,x_m];
            fx=[fx,f(x_m)];
            err=abs(f(x_m));
            xdif=[xdif,delta_x]; 
            it=it+1;
        end
end



% % Print on the screen
% if (it<nmax)
%     fprintf(' Convergence at the iteration k : %d \n',it);
% else
%     fprintf(' The maximum number of iterations was reached : %d \n',it);
% end
% fprintf(' Value calculated     : %-12.8f \n',xvect(it));

return