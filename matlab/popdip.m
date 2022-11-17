function [xk,lamk,xklist,lamklist,muklist] = popdip(x0,f,tol,maxiters,theta,kappabar,samealpha)
% POPDIP  POsitive-variables Primal-Dual Interior Point method.  This is a
% version of Algorithm 16.1 in section 16.7 of Griva, Nash, Sofer (2009), but
% specialized to
%   min         f(x)
%   subject to  x >= 0
% Uses mu and kappa formulas in section 16.7.2.  It appears the convergence
% of this algorithm it quadratic; it should be governed by Theorem 16.17.
% This implementation does not use back-tracking; only the ratio test
% for positivity, both primal and dual, adjusts the step size.  This
% implementation does not exploit any sparsity; the indefinite Newton step
% equations are solved by O(n^3) Gauss elimination.
% See documentation doc.pdf in doc/.
% Basic usage:
%   >> [xk,lamk,xklist,lamklist,muklist] = popdip(x0,f,tol,maxiters)
% where
%   x0        initial iterate
%   f         function with signature "[f,df,Hf] = f(x)" producing
%             f(x), gradient df(x), and Hessian Hf(x)
%   tol       convergence tolerance (default 10^-4)
%   maxiters  stop after this many iterations
% Examples: SMALL and OBSTACLE

    if nargin < 3,  tol = 1.0e-4;    end
    if nargin < 4,  maxiters = 200;  end
    if nargin < 5,  theta = 0.1;     end
    if nargin < 6,  kappabar = 0.9;  end
    if nargin < 7,  samealpha = true;  end

    if any(x0 <= 0.0),  error('initial iterate must be strictly feasible'),  end

    xk = x0(:);
    n = length(xk);

    % initialize dual variables:  if some components satisfy  (grad f(x0))_i > 0
    % then average to generate a scale  mu0;  otherwise guess  mu0 = 1
    [tmp, y0] = f(xk);  y0 = y0(:);  % force into column shape
    if length(y0) ~= n,  error('gradient f(x) must be same size as x'),  end
    if all(y0 <= 0)
        mu0 = 1.0;
    else
        mu0 = sum(y0(y0 > 0) .* xk(y0 > 0)) / n;  % on average:  y0*x0 = mu0
    end
    lamk = mu0 ./ xk;

    % start output lists if requested
    if nargout > 2
        xklist = [xk];
        lamklist = [lamk];
        muklist = [];
    end

    % loop: primal-dual interior point method uses Newton steps on barrier equations
    for k = 1:maxiters
        if any(xk <= 0.0)
            error('primal feasibility violated at iteration %d',k)
        end
        if any(lamk <= 0.0)
            error('dual feasibility violated at iteration %d',k)
        end
        [fxk, dfxk, Hfxk] = f(xk);          % call user's function
        meritk = merit(xk,lamk,dfxk);
        if k == 1
            merit0 = meritk;
        else
            if meritk/merit0 < tol
                break
            end
        end
        mu = min(theta*meritk,meritk^2);    % formula page 646
        M = [Hfxk,       -eye(n,n);
             diag(lamk), diag(xk)];
        %cond(M)
        c = [-dfxk + lamk;
             - lamk .* xk + mu];
        p = M \ c;                          % Gaussian elimination: O(n^3)
        kappa = max(kappabar,1.0-meritk);   % formula page 646
        if samealpha
            alpha = ratiotest([xk;lamk],p,kappa);
            xk = xk + alpha * p(1:n);
            lamk = lamk + alpha * p(n+1:2*n);
        else
            xk = xk + ratiotest(xk,p(1:n),kappa) * p(1:n);
            lamk = lamk + ratiotest(lamk,p(n+1:2*n),kappa) * p(n+1:2*n);
        end
        if nargout > 2
            xklist = [xklist xk];
            lamklist = [lamklist lamk];
            muklist = [muklist mu];
        end
    end
end

    function z = merit(x,lam,dfx)
    % MERIT implements a formula on page 642
        z = max(norm(dfx-lam),norm(lam.*x));
    end

    function alpha = ratiotest(x,dx,kappa)
    % RATIOTEST implements a formula on page 642
        alpha = 1.0;
        for j = 1:length(x)
            if dx(j) < 0
                alpha = min(alpha, - kappa * x(j) / dx(j));
            end
        end
    end
