function ivalues = polyharmonic_splines(M, L, centers, values, n, ToInterp)
% PolyharmonicSpline: Polyharmonic splines for K-dimensional scattered data (like "3D spline" in Cartool software)
%
% Polyharmonic splines are radial basis functions used for function approximation and data interpolation.
% They are very useful for interpolation of scattered data in many dimensions.
% A special case (K = 2) are thin plate splines (see tpaps.m).
%
% [Interpolated,S] = polyharmonic_splines(K, centers, values, s, centers2interp)
%
%  Inputs
% --------
% M: basis function of polyharmonic splines, based on pairwise euclidean
%    distances between centers (square form), externalized in
%    dual_interp_3D_lin_bad_elec.m to save computational time
%    externalized to save computational time
% L: pseudo-inverse of matrix containing coordinates of centers used for
%    building polyharmonic splines, ones and zeros, see
%    dual_interp_3D_lin_bad_elec.m
% centers: coordinates of values to interpolate
% values: signal fitted externally, see dual_interp_3D_lin_bad_elec.m
% n: number of dimensions in space (2 for 2D, 3 for 3D)
% ToInterp: indices of coordinates whose values should be interpolated in
%    the entire space
%
%  Outputs
% ---------
% ivalues: interpolated values
%
%-------------------------------------------------------------------------
% The MIT License (MIT)
% Copyright (c) 2015 Luke Stagner
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%-------------------------------------------------------------------------
% NB: this function is EXPERIMENTAL!
% 2019-05-09: rewritten for code optimization:
%             - removed redundant part of function (subfunction was called
%             3 times!)
%             - uses now pdist and applies operations on matrices instead
%             of embedded for loops calling subfunction (approximately 5 times faster)
%             - ~ 50% faster than the corresponding MEX file that does not
%             use pdist, squareform and matrix operations (polyharmonic_splines_4C_mex)
% 2019-05-13: full rewriting, with computation of polyharmonic splines
%             externalized as much as possible, including pdist with
%             squareform, and pinv part. polyharmonic_splines.m now only
%             interpolates signal for bad channels based on previously
%             computed weights matrix. Distance matrix is calculated only
%             once, inversion operates on good channels only (otherwise
%             the output is equivalent to the input), and
%             polyharmonic_splines.m only extract corresponding distances
%             between bad and good channels from the full matrix and multiplies
%             with previously estimated weights. Computational time has
%             been drastically reduced, speed up is very large: less than a
%             minute instead of ~ 1 hour (several days initially!)!
%-------------------------------------------------------------------------
% code refactored in Matlab language by Renaud Marquis, FBMlab, September 2018
% original Julia code found on https://gist.github.com/lstagner//04a05b120e0be7de9915#file-polyharmonic_splines-jl
%-------------------------------------------------------------------------

% if nargin < 5
%     s = 0;
% end
% [m,n] = size(centers);
% if m ~= length(values)
%     error('Mismatch between input dimensions of centers and values')
% end

%--------------- #RM@FBMlab: faster alternative ---------------
% First input argument to polyharmonicK subfunction "norm(centers(i,:) -
% centers(j,:))" for each pair of element is just the pairwise distance, so
% pdist should be more optimal here. This avoids using loops - which are
% slow in Matlab - and we can then apply other operations on the whole
% array simultaneously: this should be more efficient.
% M = squareform(pdist(XYZ)); % this is slightly different from norm(centers(i,:) - centers(j,:)), why?
% if mod(K,2)
%     M(M>=1) = M(M>=1).^K.*log(M(M>=1)); % (r^K)*log(r);
%     M(M<1) = (M(M<1).^(K-1)).*log(M(M<1).^M(M<1)); % (r.^(K-1))*log(r.^r);
% else
%     M = M.^K;
% end
% M(logical(eye(m)))=0;
% => this is now done outside of this function

% Curiously, the matrix M as obtained from above and from the original
% version (see commented lines below) are slightly different. When doing
% the calculation using sqrt( (x11-x12)^2 + (x21-x22)^2 + (x31-x32)^2 ),
% this matches the output of pdist. But when using norm, the results are
% slightly different compared to both other methods.

% IMPORTANT NOTE: Althought the difference is very small (between -5 and 5
% x 10^14), this has important consequences later, when multiplying
% coefficients (inverse "w") by M (order of magnitude in the range of -0.2
% to -2.2 in my tests, for a signal between -7 and +21 (thus ~ 7% of the
% range of the signal...)

% This occurs on R2016a at least... I assume here that results of pdist and
% manual calculations are correct. In my tests, all other variables are strictly
% equal, even "w" and "L", which are derived from "M".
% N = [ones(m,1),centers];
% 
% %--------------------------------------------------------------
% 
% % % #RM@FBMlab: the original version is slower, commenting it
% % M = zeros(m,m);
% % N = zeros(m,n+1);
% % 
% % for i = 1:m
% %     N(i,1) = 1;
% %     N(i,2:end) = centers(i,:);
% %     for j = 1:m
% %         M(i,j) = polyharmonicK(norm(centers(i,:) - centers(j,:)),K);
% %     end
% % end
% % % #RM@FBMlab: the original version is slower and commented
% 
% M = M + s.*eye(m);
% L = [[M2,N];[N' zeros(n+1,n+1)]];
% => done outside to save computational time
% => pinv(L) also done outside!
if size(centers,2)~=n
    error('Number of dimensions does not match')
end

w = L*[values;zeros(n+1,1)];
m = length(ToInterp);
ivalues = zeros(m,1);
for i = 1:m
    tmp = 0;
    l = length(w)-(n+1);
    for j = 1:l
        tmp = tmp + w(j)*M(ToInterp(i),j);
%         tmp = tmp + w(j)*Mgood(i,j);
    end
    tmp = tmp + w(l+1);
    for j = 2:n+1
        tmp = tmp + w(l+j)*centers(i,j-1);
    end
    ivalues(i) = tmp;
end

end

% % #RM@FBMlab: not needed anymore because faster alternative using pdist
% % and operations on arrays
% function Out = polyharmonicK(r,K)
% if mod(K,2)
%     if r < 1
%         Out = (r.^(K-1))*log(r.^r);
%     else
%         Out = (r^K)*log(r);
%     end
% else
%     Out = r^K;
% end
% end
