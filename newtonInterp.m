% newtonInterp takes a set of n+1 points (f1, t1), (f2, t2),..., (fn, tn)  
% and returns a polynomial, p(t), of degree at most n such that p(ti) = fi
% for i = 0, 1, 2,..., n
% tValues and fValues must be of equal length.

function [ interpolant ] = newtonInterp(tValues, fValues)
    syms x
    A = ones(length(tValues),length(tValues));
    A = tril(A);
    
    %Set up the coefficient matrix for newton interpolation
    for i = 1:length(tValues)
        for j = 1:i
            for k = 2:j
                A(i,j) = A(i,j) * (tValues(i) - tValues(k - 1));
            end
        end
    end
    
    %Use gaussian elimination to solve for the basis function weights
    basisWeights = GEWithPartialPivot(A,fValues.');

    interpolant = 0;
    
    %Contruct the interpolating polynomial from basis functions and weights
    for i = 1:length(basisWeights)
        fi = basisWeights(i);
        for j = 2:i
            fi = fi * (x - tValues(j - 1));
        end
        interpolant = interpolant + fi;
    end
            

end
