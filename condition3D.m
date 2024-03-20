%% condition3: Conditions a set of 3D points.
%
% This function conditions a set of 3D points by transforming them to a
% canonical form.

function [conditioningMatrix, transformedPoints] = condition3D(inputPoints)
    % Calculate the centroid of the input points
    centroid = mean(inputPoints(1:3, :), 2);

    % Translate points to have the centroid at the origin
    inputPoints(1:3, :) = inputPoints(1:3, :) - centroid;

    % Calculate scale such that the mean distance from the origin is sqrt(3)
    distances = sqrt(sum(inputPoints(1:3, :).^2, 1));
    scale = sqrt(3) / mean(distances);

    % Create the conditioning matrix
    conditioningMatrix = [
        scale, 0, 0, -scale * centroid(1);
        0, scale, 0, -scale * centroid(2);
        0, 0, scale, -scale * centroid(3);
        0, 0, 0, 1
    ];

    % Apply the conditioning matrix
    transformedPoints = conditioningMatrix * inputPoints;
end
