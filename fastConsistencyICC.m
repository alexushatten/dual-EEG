function ICC = fastConsistencyICC( x, y )
% fastICC: This function computes the consistency intraclass correlation
% coefficient between vectors x and y
%
% ICC = fastConsistencyICC( x, y )
%
%  Inputs
% --------
% x: n x 1 vector of observations
% y: n x 1 vector of observations
%
%  Outputs
% ---------
% ICC: consistency ICC
%
%-------------------------------------------------------------------------
% Renaud Marquis @FBMlab, June 2018
% based on Liber Eleutherios's fastcorr.c
% and Arash Salarian's ICC.m
%-------------------------------------------------------------------------

ICC = fastICC(x,y);

end

