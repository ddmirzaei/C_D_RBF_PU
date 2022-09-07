function [ps,pst] = LagrangMat(Y,X,Xt,xc,rho,op,ScaleOrd)
% This function computes the Lagrange functions of Hermite interpolation 
%  for polyharmonic splines (PHS) kernels by applying the scaling rule 
% Inputs:
%   Y: evaluation points 
%   X: trial points (centers)
%   Xt: Hermite points
%   xc: the center of stencil
%   rho: the size of stencil
%   op: operator
%   ScaleOrd: the scaling order of op
% Outputs:
%   ps: Lagrange matrix $\psi$ at Xe bases on X
%   pst: Lagrange matrix $\tild\psi$ at Xe based on XH

% points are shifted and scaled
X = (X-xc)/rho;
Xt = (Xt-xc)/rho;
Y = (Y-xc)/rho;

A = KerMat(X,X, '1');
AL = KerMat(X,Xt, op);
P = PolyMat(X,'1');
np = size(P,2); nx = size(X,1);
switch op
    case('1')
        R = KerMat(Y,X,op);
        S = PolyMat(Y,op);
        U = [A P;P' zeros(np)]\[R';S'];
        
        % re-scaling
        ps = U(1:nx,:)/(rho^ScaleOrd);
        pst = 0; 
        return;
    case ('L')
        ALL = KerMat(Xt,Xt,'L2');
        PL = PolyMat(Xt,'L');
        RL = KerMat(Y,X,'L');
        RLL = KerMat(Y,Xt,'L2');
        SL = PolyMat(Y,'L');
    case ('x')
        ALL = KerMat(Xt,Xt,'xx');
        PL = PolyMat(Xt,'x');
        RL = KerMat(Y,X,'x');
        RLL = KerMat(Y,Xt,'xx');
        SL = PolyMat(Y,'x');
    case ('y')
        ALL = KerMat(Xt,Xt,'yy');
        PL = PolyMat(Xt,'y');
        RL = KerMat(Y,X,'y');
        RLL = KerMat(Y,Xt,'yy');
        SL = PolyMat(Y,'y');
end
U = [A AL P;AL' ALL PL;P' PL' zeros(np)]\[RL';RLL';SL'];

% re-scaling
ps = U(1:nx,:)/(rho^ScaleOrd);
pst = U(nx+1:end-np,:)/(rho^(ScaleOrd-2));
