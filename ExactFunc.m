function u = ExactFunc(X,op)
% This function computes the Franke's function, its Laplacian and
%   its first derivatives at set point X
% Inputs:
%    X: test points
%    op: operator
% Outputs:
% u: vector values
global ExampleNo
x = X(:,1); y = X(:,2);
switch ExampleNo
    case 'ex1'   % Franke's function
        switch op
            case('1') % the function itself
                s1 = 3/4*exp(-1/4*((9*x-2).^2+(9*y-2).^2));
                s2 = 3/4*exp(-1/49*(9*x+1).^2-1/10*(9*y+1).^2);
                s3 = 1/2*exp(-1/4*((9*x-7).^2+(9*y-3).^2));
                s4 = 1/5*exp(-(9*x-4).^2-(9*y-7).^2);
                u = s1+s2+s3-s4;
            case('x') % its first derivative with respect to x
                s1x = 3/4*(-9/2)*(9*x-2).*exp(-1/4*((9*x-2).^2+(9*y-2).^2));
                s2x = 3/4*(-18/49)*(9*x+1).*exp(-1/49*(9*x+1).^2-1/10*(9*y+1).^2);
                s3x = 1/2*(-9/2)*(9*x-7).*exp(-1/4*((9*x-7).^2+(9*y-3).^2));
                s4x = 1/5*(-18)*(9*x-4).*exp(-(9*x-4).^2-(9*y-7).^2);
                u = s1x+s2x+s3x-s4x;
            case('y') % its first derivative with respect to y
                s1y = 3/4*(-9/2)*(9*y-2).*exp(-1/4*((9*x-2).^2+(9*y-2).^2));
                s2y = 3/4*(-9/5)*(9*y+1).*exp(-1/49*(9*x+1).^2-1/10*(9*y+1).^2);
                s3y = 1/2*(-9/2)*(9*y-3).*exp(-1/4*((9*x-7).^2+(9*y-3).^2));
                s4y = 1/5*(-18)*(9*y-7).*exp(-(9*x-4).^2-(9*y-7).^2);
                u = s1y+s2y+s3y-s4y;
            case('L') % its Laplacian
                u = 324./(5*exp((9*x - 4).^2 + (9*y - 7).^2)) ...
                    - 243./(4*exp((9*x - 2).^2/4 + (9*y - 2).^2/4))...
                    - 81./(2*exp((9*x - 7).^2/4 + (9*y - 3).^2/4))...
                    - 14337./(980*exp((9*x + 1).^2/49+(9*y + 1).^2/10))...
                    + (3*((81*x)/2 - 9).^2)./(4*exp((9*x - 2).^2/4 ...
                    + (9*y - 2).^2/4))+ ((81*x)/2 - 63/2).^2./(2*exp((9*x - 7).^2/4 ...
                    + (9*y - 3).^2/4)) - (162*x - 72).^2./(5*exp((9*x - 4).^2 ...
                    + (9*y - 7).^2)) + (3*((162*x)./49 ...
                    + 18/49).^2)./(4*exp((9*x + 1).^2/49 + (9*y + 1).^2/10)) ...
                    + (3*((81*y)/2 - 9).^2)./(4*exp((9*x - 2).^2/4 + (9*y - 2).^2/4)) ...
                    + ((81*y)/2 - 27/2).^2./(2*exp((9*x - 7).^2/4 + (9*y - 3).^2/4)) ...
                    + (3*((81*y)/5 + 9/5).^2)./(4*exp((9*x + 1).^2/49 + (9*y + 1).^2/10))...
                    - (162*y - 126).^2./(5*exp((9*x - 4).^2 + (9*y - 7).^2));                
        end
    case 'ex2'  % a cosine-form example
        switch op
            case('1') % the function itself
                u=cos(pi*x).*cos(pi*y); 
            case('x') % its first derivative with respect to x
                u= -pi*cos(pi*y).*sin(pi*x);
            case('y') % its first derivative with respect to y
                u= -pi*sin(pi*y).*cos(pi*x);
            case('L') % its Laplacian
                u = -2*pi^2*cos(pi*y).*cos(pi*x);
        end
end