function solutionVector = computeLeastSquaresSolution(A)
    % Compute the SVD of the matrix A
    [~, ~, V] = svd(A);

    % The solution is the last column of V
    solutionVector = V(:, end);

    % Normalize the solution vector if it's not already normalized
    if solutionVector(end) ~= 1
        solutionVector = solutionVector / solutionVector(end);
    end
end
