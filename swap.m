function [A, B] = swap( A, B )
% swap: stupid function that swaps variables A and B, just to avoid
% multiple lines of code...
%
% [A, B] = swap( A, B )
%
%  Inputs
% --------
% A: any variable
% B: any variable
%
%  Outputs
% ---------
% A: B
% B: A
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

Abkp = A;
A = B;
B = Abkp;

end

