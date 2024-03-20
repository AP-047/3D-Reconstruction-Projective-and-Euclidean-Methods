%% main_reconstruction: Main file for intiating the 3D reconstruction process.
%
% This function reads image points, computes the fundamental matrix,
% defines camera matrices, triangulates points, visualizes the projective
% reconstruction, reads control points, homogenizes them, triangulates the
% control points, computes the homography matrix, applies the homography
% to get the Euclidean reconstruction, and visualizes the Euclidean
% reconstruction.
function main_reconstruction()
    % Read image points from 'bh.dat'
    fileHandle = fopen('bh.dat', 'r');
    imagePoints = fscanf(fileHandle, '%f%f%f%f', [4 inf]);  % Read four floating-point numbers per line
    fclose(fileHandle);
    
    % Extract image points into x1 and x2
    x1 = imagePoints(1:2, :);
    x2 = imagePoints(3:4, :);

    % Display the size of x1 and x2
    disp('Size of x1:');
    disp(size(x1));
    disp('Size of x2:');
    disp(size(x2));

    % Compute the fundamental matrix
    F = computeFundamentalMatrix(x1, x2);

    % Display the size of F
    disp('Size of F:');
    disp(size(F));
    disp(F);

    % Define camera matrices
    [cameraMatrix1, cameraMatrix2] = defineCameraMatrices(F);
    disp('Size of Camera Matrix 1:');
    disp(size(cameraMatrix1));
    disp('Size of Camera Matrix 2:');
    disp(size(cameraMatrix2));
    disp(cameraMatrix1);
    disp(cameraMatrix2);

    % Triangulate points
    X_projective = triangulatePoints(cameraMatrix1, cameraMatrix2, x1, x2);

    % Display the size of X_projective
    disp('Size of Projective Reconstruction (Xp):');
    disp(size(X_projective));

    % Visualize projective reconstruction
    visualize3D(X_projective, 'Projective Reconstruction');

    % Read control points from 'pp.dat'
    [x1_control, x2_control, X_euclidean] = readControlPoints('pp.dat');

    % Display the size of x1_control, x2_control, and X_euclidean
    disp('Size of x1_control:');
    disp(size(x1_control));
    disp('Size of x2_control:');
    disp(size(x2_control));
    disp('Size of X_euclidean:');
    disp(size(X_euclidean));

    % Homogenize the control points
    x1_control(3, :) = 1;
    x2_control(3, :) = 1;
    X_euclidean(4, :) = 1;

    % Display the size of homogenized variables
    disp('Size of Homogenized x1_control:');
    disp(size(x1_control));
    disp('Size of Homogenized x2_control:');
    disp(size(x2_control));
    disp('Size of Homogenized X_euclidean:');
    disp(size(X_euclidean));

    % Triangulate control points
    X_euclidean_estimated = triangulatePoints(cameraMatrix1, cameraMatrix2, x1_control, x2_control);

    % Display the size of X_euclidean_estimated
    disp('Size of Euclidean Reconstruction (Xe_estimated):');
    disp(size(X_euclidean_estimated));

    % Compute the homography matrix
    homographyMatrix = computeHomography(X_euclidean_estimated, X_euclidean);

    % Display the size of homographyMatrix
    disp('Size of Homography Matrix (H):');
    disp(size(homographyMatrix));

    % Apply the homography to get the Euclidean reconstruction
    X_euclidean_reconstructed = ennorm(homographyMatrix * X_projective);

    % Display the size of X_euclidean_reconstructed
    disp('Size of Reconstructed Euclidean Points (Xe_reconstructed):');
    disp(size(X_euclidean_reconstructed));

    % Visualize Euclidean reconstruction
    visualize3D(X_euclidean_reconstructed, 'Euclidean Reconstruction');
end
%% 

% computeFundamentalMatrix: Computes the fundamental matrix from image points.
%
% This function takes two sets of corresponding image points and calculates
% the fundamental matrix using Direct Linear Transformation (DLT).

function fundamentalMatrix = computeFundamentalMatrix(imagePoints1, imagePoints2)
    % Homogenize the image points
    imagePoints1(3, :) = 1;
    imagePoints2(3, :) = 1;
    
    % Normalize the points
    [normalizedPoints1, normalizationMatrix1] = normalizePoints(imagePoints1);
    [normalizedPoints2, normalizationMatrix2] = normalizePoints(imagePoints2);
    
    % Display the sizes of normalizedPoints1 and normalizedPoints2
    disp('Size of normalizedPoints1:');
    disp(size(normalizedPoints1));
    disp('Size of normalizedPoints2:');
    disp(size(normalizedPoints2));

    % Build the matrix A for the linear system
    A = [];
    for i = 1:size(normalizedPoints1, 2)
        A = [A; normalizedPoints2(1, i)*normalizedPoints1(:, i)' normalizedPoints2(2, i)*normalizedPoints1(:, i)' normalizedPoints2(3, i)*normalizedPoints1(:, i)'];
    end
    disp('Matrix A:');
    disp(A);

    % Solve the linear system using Direct Linear Transformation (DLT)
    f = computeLeastSquaresSolution(A);
    
    % Reshape the solution into a 3x3 matrix
    fundamentalMatrix = reshape(f, 3, 3)';
    
    % Enforce rank-2 constraint on the fundamental matrix
    [U, D, V] = svd(fundamentalMatrix);
    D(3, 3) = 0;
    fundamentalMatrix = U * D * V';
    
    % Denormalize the fundamental matrix
    fundamentalMatrix = normalizationMatrix2' * fundamentalMatrix * normalizationMatrix1;
end

%% % defineCameraMatrices: Defines camera matrices from the fundamental matrix.
%
% Given the fundamental matrix, this function computes the epipoles and
% defines the camera matrices for the first and second images.

function [cameraMatrix1, cameraMatrix2] = defineCameraMatrices(fundamentalMatrix)
    % Define the first camera matrix in canonical form
    cameraMatrix1 = [eye(3), zeros(3, 1)];

    % Compute the right epipole (e2) from the fundamental matrix
    [U, D, V] = svd(fundamentalMatrix);
    epipole1 = V(:, 3);
    epipole2 = U(:, 3); % Ensure the epipole is in homogeneous coordinates

    % Display the sizes of epipole1 and epipole2
    disp('Size of epipole1:');
    disp(size(epipole1));
    disp('Size of epipole2:');
    disp(size(epipole2));

    % Define the second camera matrix
    epipole2_cross = [0, -epipole2(3), epipole2(2);
                      epipole2(3), 0, -epipole2(1);
                      -epipole2(2), epipole2(1), 0];
    cameraMatrix2 = [epipole2_cross * fundamentalMatrix + [epipole2, epipole2, epipole2], epipole2];
end

%% % triangulatePoints: Triangulates 3D points from corresponding image points.
%
% This function takes image points and camera matrices and performs
% triangulation to estimate the 3D coordinates of the points.

function triangulatedPoints = triangulatePoints(cameraMatrix1, cameraMatrix2, imagePoints1, imagePoints2)
    % Triangulates 3D points from corresponding image points and camera matrices.

    % Number of image point pairs
    numPoints = size(imagePoints1, 2);

    % Initialize the matrix to store triangulated points
    triangulatedPoints = zeros(4, numPoints);

    % Iterate over each image point pair
    for pointIndex = 1:numPoints
        % Formulate the linear system A for triangulation
        A = [
            imagePoints1(1, pointIndex) * cameraMatrix1(3, :) - cameraMatrix1(1, :);
            imagePoints1(2, pointIndex) * cameraMatrix1(3, :) - cameraMatrix1(2, :);
            imagePoints2(1, pointIndex) * cameraMatrix2(3, :) - cameraMatrix2(1, :);
            imagePoints2(2, pointIndex) * cameraMatrix2(3, :) - cameraMatrix2(2, :);
        ];

        % Solve the linear system using Direct Linear Transformation (DLT)
        triangulatedPoints(:, pointIndex) = ennorm(computeLeastSquaresSolution(A));
    end
end
%% % computeHomography: Computes the homography matrix between two sets of 3D points.
%
% This function takes two sets of 3D points, conditions them, designs a
% matrix for homography estimation, and computes the homography matrix.


function homographyMatrix = computeHomography(controlPointsXc, euclideanPointsXe)
    % Condition the control and euclidean points
    normalizationMatrix1 = condition3D(controlPointsXc);
    normalizedControlPoints = normalizationMatrix1 * controlPointsXc;
    
    normalizationMatrix2 = condition3D(euclideanPointsXe);
    normalizedEuclideanPoints = normalizationMatrix2 * euclideanPointsXe;

    % Design matrix for the homography estimation
    designMatrix = design_homo3(normalizedControlPoints, normalizedEuclideanPoints);

    % Solve the linear system using Direct Linear Transformation (DLT)
    homographyVector = computeLeastSquaresSolution(designMatrix);

    % Reshape the solution into a 4x4 matrix
    homographyMatrix = reshape(homographyVector(1:16), 4, 4)';

    % Normalize the matrix to set the last element to 1 (homogeneous coordinates)
    homographyMatrix = homographyMatrix / homographyMatrix(4, 4);

    % Undo the conditioning to obtain the final homography matrix
    homographyMatrix = inv(normalizationMatrix2) * homographyMatrix * normalizationMatrix1;
end
%% % ennorm: Homogenizes a set of 3D points.
%
% This function homogenizes a set of 3D points by dividing each coordinate
% by the last coordinate.

function y = ennorm(x)
    % Direkte Durchführung einer elementweisen Division ohne Schleife
    y = x ./ x(end, :);
end

% Weitere Hilfsfunktionen (normalizePoints, computeLeastSquaresSolution, condition3, design_homo3) wie benötigt

% visualize3D: Visualizes 3D points in a scatter plot.

% This function creates a 3D scatter plot to visualize a set of 3D points.

% Visualisierung und andere Hilfsfunktionen unten hinzufügen
function visualize3D(X, titleText, varargin)
    % Visualizes 3D points using a scatter plot with a gradient of colors.

    % Default values for optional parameters
    markerSize = 10;
    squareAxis = true;
    viewAngles = [32, 75];
    displayGrid = false;
    colorMap = parula;  % Default color map

    % Check for optional parameters
    if ~isempty(varargin)
        for i = 1:2:length(varargin)
            param = lower(varargin{i});
            value = varargin{i+1};

            % Update parameters based on input
            switch param
                case 'markersize'
                    markerSize = value;
                case 'squareaxis'
                    squareAxis = value;
                case 'viewangles'
                    viewAngles = value;
                case 'displaygrid'
                    displayGrid = value;
                case 'colormap'
                    colorMap = value;
                otherwise
                    warning('Unknown parameter: %s', param);
            end
        end
    end

    % Create a new figure
    figure;

    % Scatter plot of 3D points with color gradient
    scatter3(X(1, :), X(2, :), X(3, :), markerSize, X(3, :), 'filled', 'MarkerFaceColor', 'flat', 'CData', X(3, :));

    % Customize plot appearance
    if squareAxis
        axis square;
    end

    % Set the view angles
    view(viewAngles);

    % Set title and labels
    title(titleText);
    xlabel('X-Axis');
    ylabel('Y-Axis');
    zlabel('Z-Axis');

    % Display grid lines if specified
    if displayGrid
        grid on;
    end

    % Enhance the plot by adding a legend, if applicable
    if numel(varargin) >= 2 && strcmpi(varargin{1}, 'legend')
        legend(varargin{2:end});
    end

    % Add color bar to indicate the Z-axis values
    colorbar;
    colormap(colorMap);
    caxis([min(X(3, :)), max(X(3, :))]);  % Adjust the color range based on Z-axis values
end