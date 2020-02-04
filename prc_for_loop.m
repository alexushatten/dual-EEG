function prc_for_loop( n, N, Dot, Start, Step )
% prc_for_loop : monitor progression of for loop, 
%
%   Example usage
% -----------------
%
% for n = 1:N
%   prc_for_loop( n, N, Dot, Start, Step );
%   % processing
%   Results(n) = compute_stuff(n);
% end
%
%  Inputs
% --------
% n    : current iteration
% N    : total number of iterations
% Dot  : actual number of iterations before printing a dot ("."), default = 200
% Start: first iteration (cosmetic, to display "0%")
% Step : actual number of iterations before printing a new line, default =
%        50 * Dot, depends on the size of your Matlab command window size.
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2018
%-------------------------------------------------------------------------

if nargin < 3
    Dot = 200;
end
if nargin < 4
    Start = 1;
end
if nargin < 5
    Step = 50*Dot;
end
if rem(n,Dot)==0 && n~=Start
    fprintf('.');
end
if rem(n,Step)==0 || n==Start
    if round(n/N*100)<10
        fprintf('\n  %d%%',round(n/N*100));
    else
        fprintf('\n %d%%',round(n/N*100));
    end
end
if n == N
    fprintf('\n');
end

end

