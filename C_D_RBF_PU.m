function [C,Ct] = C_D_RBF_PU(Y,X,Xc,rho,op)
% This function produces the C and Ct matrices in 
%        (Lu)|Y = Cu|X + Ct*(Lu)|X 
%   where L is the operator using the D-RBF-PU method 
% Inputs:
%   Y: test points
%   X: trial points
%   Xc: patch centers
%   rho: patch sizes
%   op: operator
% Outputs:
%   C: the standard differentiation matrix 
%   Ct: the Hermite differentiation matrix

global WeightType  % PU weight type: Smooth or Constant-generated
global cH          % to define Hermite points in each stencil 
N = size(X,1); Ne = size(Y,1); Nc = size(Xc,1);  % sizes ...
C = spalloc(Ne,N,0);                % initial allocation
Ct = spalloc(Ne,Ne,0);              % initial allocation
ScaleOrd = ScalingOrder (op);       % scaling order of the operator
IndX = PointsInPatch(X,Xc,rho);     % indices of trial points in each patch
IndY = PointsInPatch(Y,Xc,rho);     % indices of test points in each patch
IndXc = PointsInPatch(Xc,Xc,2*rho); % indices of patch centers in each patch
switch WeightType
    case 'Smooth'
        for j=1:Nc
            % xc: patch center number j
            xc = Xc(j,:);
            
            % Xj: all trial points in patch j
            ind = IndX{j}; inde = IndY{j}; indc = IndXc{j};
            Xj = X(ind,:); Nind = length(ind);
            
            % Yj: all test points in patch j
            Yj = Y(inde,:);

            % XHj: Hermite points in patch j
            indb = find(DistMat(Yj,xc) > cH*rho);
            indH = inde(indb); NindH = length(indH);
            XHj = Y(indH,:);
            
            if isempty(inde)==0
                
                % w: PU weights at test points 
                w = PUweight(Yj,Xc(indc,:),xc,rho);
                
                % c, ct: weight matrices (Lagrance functions at Yj)
                [c,ct] = LagrangMat(Yj,Xj,XHj,xc,rho,op,ScaleOrd);
                
                % inserting weights into global matrices and updating
                C(inde,ind) = C(inde,ind) + repmat(w,1,Nind).*c';
                Ct(inde,indH) = Ct(inde,indH) + repmat(w,1,NindH).*ct';
            end
        end
    case 'ConstGen'
        for j=1:Nc
            
            % xc: patch center number j
            xc = Xc(j,:);    
            
            % Xj: all trial points in patch j
            indP = IndX{j}; inde0 = IndY{j}; indc = IndXc{j};
            Xj = X(indP,:);
                        
            % finding indices of test points in which xc is their closest center
            Yj0 = Y(inde0,:);
            D = DistMat(Yj0,Xc(indc,:));
            [Dmin,Dind] = min(D,[],2); % find minimum distances           
            inde = (Dind == 1);
            indx = inde0(inde); % going from local to global indices
            
            % XHj: Hermite points in patch j
            indb = find(DistMat(Yj0,xc) > cH*rho);
            indH = inde0(indb);
            XHj = Y(indH,:);

            if ~isempty(inde)
                % test points in tile j
                Yj = Yj0(inde,:);
                
                % c, ct: weight matrices (Lagrance functions at Yj)
                [c,ct] = LagrangMat(Yj,Xj,XHj,xc,rho,op,ScaleOrd);
                
                % inserting weights into global matrices 
                C(indx,indP) = c';
                Ct(indx,indH) = ct';
            end
        end
end
