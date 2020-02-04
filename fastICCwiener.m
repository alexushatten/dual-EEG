function icc = fastICCwiener( x, y )
% fastICCwiener: consistency ICC, aka Winer's adjustment for anchor points
%
% icc = fastICCwiener( x, y )
%
%  Inputs
% --------
% x : n x 1 vector of observation for variable X
%
% y : n x 1 vector of observation for variable Y
%
%  Outputs
% ---------
% icc : degree of consistency among measurements, also known as
%       norm-referenced reliability and as Winer's adjustment for anchor
%       points
%
%-------------------------------------------------------------------------
% based on Arash Salarian's ICC.m & Liber Eleutherios's fastcorrelation.m
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2018
%-------------------------------------------------------------------------

if exist('fastICC','file')==3
    icc = fastICC(x,y);
else
    M = [x,y];
    [n, k] = size(M);
    SStotal = var(M(:)) *(n*k - 1);
    MSR = var(mean(M, 2)) * k;
    MSC = var(mean(M, 1)) * n;
    MSE = (SStotal - MSR *(n - 1) - MSC * (k -1))/ ((n - 1) * (k - 1));
    icc = (MSR - MSE) / (MSR + (k-1)*MSE);
end

end

