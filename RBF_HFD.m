function [C,Ct] = RBF_HFD(Y,X,~,delta,op)
% This function produces the C and Ct matrices in 
%        (Lu)|Y = Cu|X + Ct*(Lu)|X 
%   where L is the operator using the RBF-HFD method 
% Inputs:
%   Y: test points
%   X: trial points
%   delta: size of stencils
%   op: operator
% Outputs:
%   C: the standard differentiation matrix 
%   Ct: the Hermite differentiation matrix

global cH    % to define Hermite points in each stencil 
N = size(X,1); Ne = size(Y,1);      % sizes ...
C = spalloc(Ne,N,0);        % initial allocation
Ct = spalloc(Ne,Ne,0);      % initial allocation
IndX = PointsInPatch(X,Y,delta);   % indices of trial points in each patch
IndY = PointsInPatch(Y,Y,delta);   % indices of trial points in each patch
ScaleOrd = ScalingOrder(op);        % scaling order of the operator
for k = 1:Ne
    yk = Y(k,:);       % test point
    ind = IndX{k};    % indices of points in stencil k
    inde = IndY{k};    % indices of points in stencil k
    Xk = X(ind,:);     % points in stencil k
    indb = find(DistMat(Y(inde),yk) > cH*delta); 
    indH = inde(indb);  % indices of Hermite points in the stencil
    XH = Y(indH,:);    % Hermite points in stencil k
    
    % weight vectors (Lagrance functions at yk)
    [c,ct] = LagrangMat(yk,Xk,XH,yk,delta,op,ScaleOrd); 
    
    % inserting weights into global matrices 
    C(k,ind) = c;
    Ct(k,indH) = ct;
end
