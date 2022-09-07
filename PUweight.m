function w = PUweight(X,Xc,xc,del)
% This function computes the Shepard weights at X (localted at a patch with
%  radius del and center xc). But all patch with centers Xc are contributed 
% Inouts:
%   X: evaluation points
%   Xc: contributed patches
%   xc: the center of that patch containing X
%   del: patch radii 
% Outputs:
%   w: weight vector
%
r0 = DistMat(X,xc)/del; r = DistMat(X,Xc)/del;
p0 = max(1-r0,0).^4.*(4*r0+1); p = max(1-r,0).^4.*(4*r+1);
w = p0./sum(p,2);
