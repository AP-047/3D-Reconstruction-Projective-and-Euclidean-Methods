%% normalizePoints: Normalizes a set of 2D points.
%
% This function normalizes a set of 2D points using a transformation matrix.

function [normalizedPts, normalizationMatrix] = normalizePoints(inputPts)
    % Calculate the centroid of the input points
    centroid = mean(inputPts(1:2, :), 2);

    % Shift the points to have the centroid at the origin
    shiftedPts = inputPts(1:2, :) - centroid;

    % Calculate the average distance from the origin
    averageDistance = mean(sqrt(sum(shiftedPts.^2, 1)));

    % Calculate the scaling factor
    scaleFactor = sqrt(2) / averageDistance;

    % Create the normalization matrix
    translationMatrix = eye(3);
    translationMatrix(1:2, 3) = -centroid;

    scalingMatrix = diag([scaleFactor, scaleFactor, 1]);

    normalizationMatrix = scalingMatrix * translationMatrix;

    % Apply the normalization transformation to the input points
    normalizedPts = normalizationMatrix * inputPts;
end
