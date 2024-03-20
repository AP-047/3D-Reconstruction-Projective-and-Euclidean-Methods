function [x1, x2, Xe] = readControlPoints(filename)
    fh = fopen(filename, 'r');
    A = fscanf(fh, '%f%f%f%f%f%f%f', [7 inf]);
    fclose(fh);
    
    % Extract data
    x1 = A(1:2, :);
    x2 = A(3:4, :);
    Xe = A(5:7, :);
    
    % Ensure they are in 3xN format
    x1(3, :) = 1;
    x2(3, :) = 1;
    Xe(4, :) = 1; % Homogenize
end
