function [ xStar, error ] = bisectSolve(f, a, b, tol)
        syms x;
        error = [];
        
        while (abs(a - b) > (tol * 2))
            bisection = a + (b - a) / 2; % to avoid underflow error
            error = [error bisection];
            if (sign(subs(f,bisection)) == sign(subs(f,a)) )
                a = bisection;
            else
                b = bisection;
            end
        end
        
        xStar = a + (b - a) / 2;
        error = [error xStar];
        error = abs(error - xStar);
        
end