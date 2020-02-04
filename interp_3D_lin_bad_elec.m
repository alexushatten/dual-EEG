function EEG = interp_3D_lin_bad_elec( EEG, XYZ, BadChanIdx, Method )
% dual_interp_3D_lin_bad_elec: interpolates bad EEG channels as in Cartool
% using 3D linear interpolation
%
% EEG = dual_interp_3D_lin_bad_elec( EEG, XYZ, BadChanIdx, Method )
%
%  Inputs
% --------
% EEG: [channel x time] EEG traces
% XYZ : [N x 3] electrode coordinates in x, y, z axes, with normalized coordinates *
% BadChanIdx: vector of integer, indices of bad electrode to interpolate
% Method:   'polyharmonic_splines': very fast, close to Cartool's output,
%           can be extended to N-dimensions (and potentially to higher
%           spline degrees)
%           'scattered_interpolant': slow (~ 13.5x slower than polyharmonic splines),
%           very close to Cartool's output, is restricted to linear
%           (or nearest neighbour) interpolation for 2-D and 3-D
%
%  Outputs
% ---------
% EEG: EEG traces with bad channels interpolated
%
%-------------------------------------------------------------------------
% * see Cartool reference guide (Tracks interpolation - Technical points)
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
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
%             minute instead of ~ 1 hour (several days initially!)! In the
%             end, this function is 165 x faster than Cartool on a
%             quad-core PC. We can thus still expect a speed-up of at least
%             40 orders of magnitude when adjusting for the number of cores.
% 2019-05-14: outputs of Cartool vs. this code are almost identical, except
%             for an offset. Given these results, THIS FUNCTION IS STILL
%             EXPERIMENTAL, USE AT YOUR OWN RISK!
% IMPORTANT: THIS FUNCTION IS NOT APPROPRIATE WHEN CHANGING SPACE (e.g.
%            256 -> 128 channels, only for interpolating electrodes while
%            staying in the same electrode coordinates space)
%-------------------------------------------------------------------------
% code refactored in Matlab language by Renaud Marquis, FBMlab, September 2018
% original Julia code found on https://gist.github.com/lstagner//04a05b120e0be7de9915#file-polyharmonic_splines-jl
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, October 2018, last updated May 2019
%-------------------------------------------------------------------------
% The MIT License (MIT)
% Copyright (c) 2015 Luke Stagner
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%-------------------------------------------------------------------------

tic;
GoodChanIdx = setdiff((1:size(EEG,1))',BadChanIdx);

switch lower(Method)
    case 'polyharmonic_splines'
        
        % Although with m >=3 is sufficient according to [1],
        % here we will use a spline degree of 1 (linear), because it will
        % ensure that the data the interpolated channels will never exceed the data
        % minimuma and maximum, which is important to avoid bad surprises if
        % fitting goes wrong...
        % [1] Ferree, T. C. (2000). Spline Interpolation of the Scalp EEG., www.egi.com/Technotes/SplineInterpolation.pdf
        % The degree of the spline is equal to (m - 1), which means that with m = 2
        % we get linear interpolation:
        K = 2;
        S = 0; % smoothing (off by default)
        XYZ_good_chans = XYZ(GoodChanIdx,:);
        
        % GoodChanIdx = setdiff(1:size(EEG,1),BadChanIdx);
        % XYZ = XYZ(GoodChanIdx,:);
        [m,n] = size(XYZ);
        % we calculate distance matrix only once because it does not change:
        M = squareform(pdist(XYZ)); % this is slightly different from norm(centers(i,:) - centers(j,:)), why?
        % and we also make the polyharmonic spline
        if mod(K,2)
            M(M>=1) = M(M>=1).^K.*log(M(M>=1)); % (r^K)*log(r);
            M(M<1) = (M(M<1).^(K-1)).*log(M(M<1).^M(M<1)); % (r.^(K-1))*log(r.^r);
        else
            M = M.^K;
        end
        M(logical(eye(m)))=0;
        Mgood = M(GoodChanIdx,GoodChanIdx); M = M(:,GoodChanIdx);
        N = [ones(length(GoodChanIdx),1),XYZ_good_chans];
        Mgood = Mgood + S.*eye(length(GoodChanIdx)); % smoothing would occur here
        L = [[Mgood,N];[N' zeros(n+1,n+1)]];
        Linv = pinv(L);
        % doing all time points with the estimated distance matrix:
        for t = 1:size(EEG,2)
            Dinterp = polyharmonic_splines(M,Linv,XYZ(BadChanIdx,:),EEG(GoodChanIdx,t),n,BadChanIdx);
            % still 50% faster than MEX version without pdist & squareform (unsupported) and less matrix operations
            % 2019-05-12: even after having externalized pdist and squareform and
            % made only the rest in MEX form, the non-MEX is still faster...
            %     Dinterp = polyharmonic_splines_mex(3,M,XYZ,EEG(:,t),1); % still 50% faster than MEX version without pdist & squareform (unsupported) and less matrix operations
            EEG(BadChanIdx,t) = Dinterp;
        end
        
    case 'scattered_interpolant'
        %% Alternative method
        % - works only for 2D and 3D, and only with linear interpolation
        % - slower
        % - slightly closer to Cartool output
        
        % evaluate scattered interpolant only once
        F = scatteredInterpolant(XYZ(GoodChanIdx,:),EEG(GoodChanIdx,1),'linear');
        for t = 1:size(EEG,2)
            prc_for_loop(t,size(EEG,2),100);
            % just update values, more efficient
            F.Values = EEG(GoodChanIdx,t);
            % query
            Dinterp = F(XYZ(BadChanIdx,:));
            % update
            EEG(BadChanIdx,t) = Dinterp;
        end
        
    otherwise
        error('Invalid method string')
end
Toc = toc;
if Toc<60
    fprintf('3D spline (linear) interpolation\nperformed in %d seconds.\n',round(Toc));
elseif (Toc/60)<60
    fprintf('3D spline (linear) interpolation\nperformed in %d minutes.\n',round(Toc/60));
else
    fprintf('3D spline (linear) interpolation\nperformed in %d hours and %d minutes.\n',round(Toc/3600),round(rem(Toc,3600)/60));
end
end

