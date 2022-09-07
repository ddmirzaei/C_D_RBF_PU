%%   
% Execut this script to see the results of paper
%  S. Arefian, D. Mirzari, A compact radial basis function partition of
%     unity method, Comput. Math. Appl. 2022. 
% We consider an equation of the form ((-kappa*Delta u + gamma*u = f ))
% with Dirichlet and Neumann boundary conditions on the unit square. 
%%
clc
clear

global kappa gamma ExampleNo
disp('------------- Results of the paper --------------------')
% Example 1 here is the example given in the above paper, 
% the Poisson equation with Franke's function as an exact solution
ExampleNo = 'ex1'; 
kappa = 1; gamma = 0;
disp('-------------------------------------------------------')
disp('Compact D-RBF-PU method strats ...')
disp('WeithType = Smooth, RBFtype = r^k log r, (k = 4,6,8), BC = Dirichlet')
WeightType = 'Smooth'; BoundaryType = 'Dirichlet'; PointType = 'halton';                                  
RBFtype = 'tp'; MethodType = 'pu'; pVec = [1 2 3]; 
[N,InfErr,Orders,SetupTime,SolveTime,nz_percent] = ...
      PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,MethodType,...
          WeightType,BoundaryType,RBFtype,pVec)

disp('-------------------------------------------------------')
disp('Compact D-RBF-PU method strats ...')
disp('WeithType = Smooth, RBFtype = r^b, (b = 5,7,9), BC = Dirichlet')
WeightType = 'Smooth'; BoundaryType = 'Dirichlet'; PointType = 'halton';                                  
RBFtype = 'p'; MethodType = 'pu'; pVec = [1 2 3];
[N,InfErr,Orders,SetupTime,SolveTime,nz_percent] = ...
      PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,MethodType,...
          WeightType,BoundaryType,RBFtype,pVec)

disp('-------------------------------------------------------')
disp('Compact D-RBF-PU method strats ...')
disp('WeithType = Constant-generated, RBFtype = r^k log r (k = 4,6,8), BC = Dirichlet')
WeightType = 'ConstGen'; BoundaryType = 'Dirichlet'; PointType = 'halton';                                  
RBFtype = 'tp'; MethodType = 'pu'; pVec = [1 2 3];  
[N,InfErr,Orders,SetupTime,SolveTime,nz_percent] = ...
      PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,MethodType,...
          WeightType,BoundaryType,RBFtype,pVec)

disp('-------------------------------------------------------')
disp('Compact D-RBF-PU method strats ...')
disp('WeithType = Constant-generated, RBFtype = r^b, (b = 5,7,9), BC = Dirichlet')
WeightType = 'ConstGen'; BoundaryType = 'Dirichlet'; PointType = 'halton';                                  
RBFtype = 'p'; MethodType = 'pu'; pVec = [1 2 3];  
[N,InfErr,Orders,SetupTime,SolveTime,nz_percent] = ...
      PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,MethodType,...
          WeightType,BoundaryType,RBFtype,pVec)

disp('-------------------------------------------------------')
disp('Compact D-RBF-PU method strats ...')
disp('WeithType = Smooth, RBFtype = r^k log r (k = 4,6,8), BC = Neumann')
WeightType = 'Smooth'; BoundaryType = 'Neumann'; PointType = 'halton';                                  
RBFtype = 'tp'; MethodType = 'pu'; pVec = [1 2 3];    
[N,InfErr,Orders,SetupTime,SolveTime,nz_percent] = ....
      PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,MethodType,...
          WeightType,BoundaryType,RBFtype,pVec)

disp('-------------------------------------------------------')
disp('Compact D-RBF-PU method strats ...')
disp('WeithType = Smooth, RBFtype = r^b, (b = 5,7,9), BC = Neumann')
WeightType = 'Smooth'; BoundaryType = 'Neumann'; PointType = 'halton';                                  
RBFtype = 'p'; MethodType = 'pu'; pVec = [1 2 3];  
[N,InfErr,Orders,SetupTime,SolveTime,nz_percent] = ...
      PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,MethodType,...
          WeightType,BoundaryType,RBFtype,pVec)

disp('-------------------------------------------------------')
disp('Comparison of RBF-HFD and Compact D-RBF-PU method strats ...')
disp('RBFtype = r^7, BC = Neumann')
WeightType = 'Smooth'; BoundaryType = 'Neumann'; PointType = 'halton';                                  
RBFtype = 'p'; MethodType = 'pu'; pVec = 2; 
[N,InfErr(:,1),Orders(1),SetupTime(:,1),SolveTime(:,1),nz_percent(:,1)] = ...
    PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
WeightType = 'ConstGen'; 
[N,InfErr(:,2),Orders(2),SetupTime(:,2),SolveTime(:,2),nz_percent(:,2)] = ...
    PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
MethodType = 'fd'; 
[N,InfErr(:,3),Orders(3),SetupTime(:,3),SolveTime(:,3),nz_percent(:,3)] = ...
    PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);

Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,MethodType,...
          WeightType,BoundaryType,RBFtype,pVec)
disp('-------------------------------------------------------')

%%
disp('------------- Additional Results --------------------')
% Example 2 is an extra example, with kappa = 1, gamma = 1/2
% 
ExampleNo = 'ex2'; 
kappa = 1/2; gamma = 1;
disp('-------------------------------------------------------')
disp('Compact D-RBF-PU method strats ...')
disp('WeithType = Smooth, RBFtype = r^k log r, (k = 4,6,8), BC = Neumann')
WeightType = 'Smooth'; BoundaryType = 'Neumann'; PointType = 'halton';                                  
RBFtype = 'tp'; MethodType = 'pu'; pVec = [1 2 3]; 
[N,InfErr,Orders,SetupTime,SolveTime,nz_percent] = ...
      PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,MethodType,...
          WeightType,BoundaryType,RBFtype,pVec)

disp('-------------------------------------------------------')
disp('Compact D-RBF-PU method strats ...')
disp('WeithType = Constant-generated, RBFtype = r^k log r (k = 4,6,8), BC = Dirichlet')
WeightType = 'ConstGen'; BoundaryType = 'Dirichlet'; PointType = 'halton';                                  
RBFtype = 'tp'; MethodType = 'pu'; pVec = [1 2 3];  
[N,InfErr,Orders,SetupTime,SolveTime,nz_percent] = ...
      PDEsolver(MethodType,WeightType,BoundaryType,PointType,RBFtype,pVec);
Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,MethodType,...
          WeightType,BoundaryType,RBFtype,pVec)

disp('-------------------------------------------------------')
