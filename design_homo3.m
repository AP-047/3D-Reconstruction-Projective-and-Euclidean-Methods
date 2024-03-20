function A = design_homo3(X1, X2)
% DESIGN_HOMO3 constructs the matrix A for homography estimation
% between two sets of 2D points X1 and X2.
A = [];
for i = 1 : size(X1, 2) % For all object points
A = [ A;
-X2(4,i)*X1(:,i)' 0 0 0 0 0 0 0 0 X2(1,i)*X1(:,i)';
0 0 0 0 -X2(4,i)*X1(:,i)' 0 0 0 0 X2(2,i)*X1(:,i)';
0 0 0 0 0 0 0 0 -X2(4,i)*X1(:,i)' X2(3,i)*X1(:,i)'];
end