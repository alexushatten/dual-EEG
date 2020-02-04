function Y = sliding_win_fun( X, Window, Fun )
% sliding_win_fun: apply function in sliding-window fashion (backward +
% forward) to signal, uses symmetric padding (padarray.m)
%
% Y = sliding_win_fun( X, Window, Fun )
%
%  Inputs
% --------
% X: [n x 1] array
% Window: integer, window size
% Fun: function handle (e.g. @nanmean)
%
%  Outputs
% ---------
% Y: [n x 1] array with sliding function "Fun" is applied to X with window
% Window
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, July 2019
%-------------------------------------------------------------------------

if numel(X)~=length(X)
    error('X must be a vector (1-dimensional array)')
end

X = X(:);
HalfWin = round(Window/2);
Xpad = padarray(X,HalfWin,'symmetric');

Windows = nan(length(X),2*HalfWin+1);
for t = 1:(length(X))
    Windows(t,:) = t:(t+HalfWin*2);
end

Y = nan(size(X));
for t = 1:length(X)
    Y(t) = Fun(Xpad(Windows(t,:)));
end

end

