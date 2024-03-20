function H = computeHomography(Xc, Xe)
    % Check if Xc and Xe are both in 3xN format
    %[rowsXc, colsXc] = size(Xc);
    %[rowsXe, colsXe] = size(Xe);
    
   % if rowsXc ~= 3 || rowsXe ~= 3 || colsXc ~= colsXe
    %%end
    
    % Initialize the design matrix
    A = [];

    % Construct the design matrix
    for i = 1:colsXc
        xc = Xc(:, i)';
        xe = Xe(:, i)';
        
        A(2*i-1 : 2*i, :) = [0, 0, 0, -xe(3)*xc, xe(2)*xc; ...
                            xe(3)*xc, 0, -xe(1)*xc, 0, -xe(3)*xc; ...
                            -xe(2)*xc, xe(1)*xc, 0, 0, -xe(2)*xc];
    end

    % Compute the singular value decomposition of A
    [~, ~, V] = svd(A);

    % Extract the homography matrix from the last column of V
    H = reshape(V(:, end), 3, 3)';
end
