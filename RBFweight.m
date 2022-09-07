function rbf = RBFweight(X,Y,RBFscale,char)
r=DistMat(X,Y);
switch (char)
    case('0')
        rbf=max(1- RBFscale*r,0).^6.*(4*(RBFscale*r).^2+18*RBFscale*r+3); 
end