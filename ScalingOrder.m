function ScaleOrd = ScalingOrder (op)
% This function returns the scaling order of operator op
% Inputs:
%  op: operator
% Outputs:
%  ScaleOrd: the scaling order
%
switch op
    case '1', ScaleOrd = 0; case 'x', ScaleOrd = 1; case 'y', ScaleOrd = 1;
    case 'xx', ScaleOrd = 2;case 'xy',ScaleOrd = 2; case 'yy',ScaleOrd = 2;
    case 'L', ScaleOrd = 2;
end
