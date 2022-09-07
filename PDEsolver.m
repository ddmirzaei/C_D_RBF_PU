function [N,InfErr,Orders,SetupTime,SolveTime,nz_percent] = ...
            PDEsolver(MethodType,weightType,BoundaryType,PointType,RBFtype,pVec)
% This function solves the PDE using both RBF-HFD and compact D-RBF-PU
%   methods. It works for the Poisson equation \Delta u  =  f with 
%   Dirichlet or Neumann boundary conditions. 
%           
% Inputs:
%   MethodType: PU or FD
%   weightType: Smooth or Constant-generated
%   BoundaryType: Dirichlet or Neumann (on bottom and top)
%   PointType: grid or Halton
%   RBFtype: powers or tps (PHS kernels)
% Outputs:
%   N: number of points (vector)
%   InfErr: the computed norm-infinity errors
%   Order: computational convergence orders 
%   SetupTime: the cpu times for setting up the final system
%   SolveTime: the cpu times for solving the final system
%
global RBFinfo cH WeightType
global gamma kappa

kp = kappa; gm = gamma;
RBFinfo.type = RBFtype;
WeightType = weightType;
switch MethodType
    case 'fd'
        CallFun = @(x,y,z,w,s) RBF_HFD(x,y,z,w,s);
    case 'pu'
        CallFun = @(x,y,z,w,s) C_D_RBF_PU(x,y,z,w,s);
end
cH = 0.75;
nlevel = 8;
for p = 1:length(pVec)
    if strcmp(RBFinfo.type, 'tp')
       RBFinfo.par = 2+2*pVec(p); RBFinfo.poly = RBFinfo.par/2+1;
    else
       RBFinfo.par = 3+2*pVec(p); RBFinfo.poly =ceil(RBFinfo.par/2);  
    end
    h = 0.1;
    for k = 1:nlevel        
        % X: trial points
        % XI: internal test points
        % XB: all boundary points
        % Xb: boundary points on the bottom side of the square
        % Xt: boundary points on the top side of the square
        % Xl: boundary points on the left side of the square
        % Xr: boundary points on the right side of the square
        hh(k) = h;
        [X,XI,XB,Xb,Xt,Xl,Xr] = PointsInSquare(0,1,0,1,h,PointType);       
        N(k)=size(X,1); NI = size(XI,1); Nb = size(Xb,1); Nt = size(Xt,1);
        
        % patch centers (grid points):
        hc = 4*h; rho = hc;
        nc = ceil(sqrt(N(k)/16)); x = linspace(2*h,1-2*h,nc);
        [x,y] = meshgrid(x,x); 
        Xc = [x(:) y(:)];
        
        tic
        cH = 0.75;
        [C_L,Ct_L]=CallFun(XI,X,Xc,rho,'L');     
        if strcmp(BoundaryType,'Neumann') % % Neumann BC
            cH = 1;
            [C_N1,Ct_N1] = CallFun(Xb,X,Xc,rho,'y');
            [C_N2,Ct_N2] = CallFun(Xt,X,Xc,rho,'y');
            [C_D1] = CallFun(Xr,X,Xc,rho,'1');
            [C_D2] = CallFun(Xl,X,Xc,rho,'1');
            A = [-kp*C_L+[gm*(speye(NI)-Ct_L) zeros(NI,N(k)-NI)];C_N1;C_N2;C_D1;C_D2];
            f = -kp*ExactFunc(XI,'L')+gm*ExactFunc(XI,'1');
            g1 = ExactFunc(Xb,'y'); g2 = ExactFunc(Xt,'y');
            R = [(speye(NI)-Ct_L)*f;(speye(Nb)-Ct_N1)*g1;(speye(Nt)-Ct_N2)*g2;...
                 ExactFunc(Xr,'1');ExactFunc(Xl,'1')];              
        else % Dirichlet BC
            [C_B] = CallFun(XB,X,Xc,rho,'1');
            A = [-kp*C_L+[gm*(speye(NI)-Ct_L) zeros(NI,N(k)-NI)]; C_B];
            f = -kp*ExactFunc(XI,'L')+gm*ExactFunc(XI,'1');
            R = [(speye(NI)-Ct_L)*f;ExactFunc(XB,'1')];
        end
        SetupTime(k,p) = toc;  % setup time
        tic
        Uap = A\R; 
        SolveTime(k,p) = toc;  % solving time
        InfErr(k,p) = norm(Uap-ExactFunc(X,'1'),inf)./norm(ExactFunc(X,'1'),inf);
        nz_percent(k,p) = nnz(A)/prod(size(A))*100;
        
        h = h/sqrt(2);  % point refinement         
    end    
    % computational orders using a linear curve fitting to errors
    a = polyfit(log(sqrt(N')),log(InfErr(:,p)),1);
    Orders(p) = -a(1);
end
