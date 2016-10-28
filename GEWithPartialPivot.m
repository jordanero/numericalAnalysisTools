% GEWithPartialPivot uses Gaussian Elimination with partial pivoting to 
% solve the system of linear equations, Ax=b, where A is a square n
% dimensional coefficient matrix, and x and b are n-vectors.
% The function takes A and b as input and returns x. An error is thrown if
% A or b is of incorrect dimensions or if A is singular. 

function [x] = GEWithPartialPivot(A, b)

    % check dimensions of A and b
    n = checkDimensions(A, b);

    % create augmented matrix
    A = [A b];

    % transform A into an upper triangular matrix
    for j = 1:(n-1)
        A = partialPivot(A, j);
        A = eliminateColumn(A, j);
        checkSingularity(A, j);
    end

    x = backwardSubstitution(A);


end

function [n] = checkDimensions(A, b)

    n = size(A);

    if n ~= size(A')
        error('Matrix A is not square.')
    end

    n = n(1);

    if size(A) ~= size(b,1)
        error('Height of vector b is different than dimension of matrix A')
    end

end

function [A] = partialPivot (A, j)

    [nRow,~] = size(A);

    % find the index, i, of the largest entry of the jth colmun below 
    % or in the jth row
    [~,i] = max(abs(A(j:nRow,j)));
    i = i + j - 1;

    % create and apply the appropriate permutation matrix to switch row j and 
    % row i
    p = 1:nRow;
    p(i) = j;
    p(j) = i;
    A = A(p,:);

end

function [A] = eliminateColumn(A, j)


    [n,~] = size(A);

    % create column j, Mjj, of the elementary multiplication matrix, Msubj.
    % Msubj * A eliminates all entries of column j below the diagonal. All
    % other columns of Msubj are equal to those of an n x n identity matrix
    % and don't need to be explicitly represented.
    Mjj = -1 * A(:,j) / A(j,j);

    % Calculate Msubj * A to eliminate non-diagonal entries of column j in
    % O(n^2)
    for i = (j+1:n)
        A(i,:) = A(i,:) + Mjj(i)*A(j,:);
end

end

function [] = checkSingularity(A, j)

    if(A(j,j) == 0)
        error('error: Martrix is singular to working precision');
    end 

end

function [x] = backwardSubstitution(A)


    [n,~] = size(A);

    for j = n: -1 : 1
        for i = 1 : (j - 1)
            % the n-vector, b, is changed, while entries in the coefficient
            % are not explicitly changed to save operations
            A(i,n+1) = A(i,n+1) - (A(i,j) / A(j,j))*A(j,n+1)
        end
        A(j, n+1) = A(j, n+1)/ A(j,j)
    end
    x = A(:,n+1) 
    
end