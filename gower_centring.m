function [G,M] = gower_centring( M )
% gower_centring: Gower centre matrix M such that rows and columns sum to 0
%
% G = gower_centring( M )
%
%  Inputs
% --------
% M : matrix
%
%  Outputs
% ---------
% G : centred matrix
%
%
%-------------------------------------------------------------------------
% Based on Legendre, 2018
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab & LREN, June 2019
%-------------------------------------------------------------------------

G = (eye(size(M))-ones(size(M))./length(M))*M*(eye(size(M))-ones(size(M))./length(M));

end

