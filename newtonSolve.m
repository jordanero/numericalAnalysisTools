% newtonSolve uses Newton's method to find a single zero for a  
% univariate equation within a given threshold for error.
% 
% newtonSolve takes f, a function, initialValue, a guess for a zero
% of f and tol, the tolerance for error with respect to the independent
% variable. newtonSolve returns the zero of f.
%
% NewtonSolve may fail under the same conditions that newton's method may 
% fail. newtonSolve throws an error if iteration results in division by 
% zero, but does not throw an error if newtons iteration is caught in loop.


function [ xStar, error ] = newtonSolve(f, initialValue, tol)
    syms x
    df = diff(f);
    
    xi = initialValue;
    xiPlus1 = newtonGuess(f, df, xi);

    while (abs(xi - xiPlus1) > tol)
        xi = xiPlus1;
        xiPlus1 = newtonGuess(f, df, xi);
    end
    
    xStar = xiPlus1;

end

% newtonGuess takes f, a function, df, its derivative and a guess
% for a zero of the function and returns a guess of a zero of f using
% the iteration step of Newton's method for solving nonlinear equations
function [xiPlus1] = newtonGuess(f, df, xi)
    syms x
    if (subs(df, x, xi) == 0)
        error('Choice of intialValue results in division by 0')
    end
    
    xiPlus1 = double(xi - (subs(f, xi) / subs(df, x, xi)));
end
    
    

