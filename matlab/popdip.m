function [xk,tauk,lamk,iteratelist,nulist,mulist] = popdip(...
                 x0,f,A,b,rtol,atol,maxiters,theta,kappabar)
% POPDIP  POsitive-variables Primal-Dual Interior Point method.
% Solves the following problems:
%   min         f(x)
%   subject to  Ax = b
%               x >= 0
% See documentation doc.pdf in doc/.  This implementation does not exploit
% sparsity; the Newton step equations are solved by O(n^3) Gauss elimination.
% Compare Algorithm 16.1 in section 16.7 of Griva, Nash, Sofer (2009).
%
% Basic usage with default parameters:
%   >> [x,tau,lam] = popdip(x0,f,A,b)
% where
%   x0        initial iterate for primal variables; in R^n
%   f         function on R^n with signature "[f,df,Hf] = f(x)" producing
%             f(x), gradient df(x), and Hessian Hf(x)
%   A         m by n matrix with full row rank
%   b         right-side vector; in R^m
%   x         primal variables at solution
%   tau       equality constraint Lagrange multipliers at solution
%   lam       positivity Lagrange multipliers (dual variables) at solution
%
% Usage on a problem with no equality constraints (min f(x) s.t. x >= 0):
%   >> [x,_,lam] = popdip(x0,f,[],[])
%
% Examples: SMALL and OBSTACLE

    if nargin < 5,  rtol = 1.0e-4;    end
    if nargin < 6,  atol = 1.0e-50;    end
    if nargin < 7,  maxiters = 200;  end
    if nargin < 8,  theta = 0.1;     end
    if nargin < 9,  kappabar = 0.9;  end

    % check sizes and force column vectors
    if any(x0 <= 0.0)
        error('primal initial iterate x0 must be strictly feasible')
    end
    xk = x0(:);
    n = length(xk);
    [m, nn] = size(A);   % "m>0" is condition for whether equality
                         % constraints are present
    if m > 0
        if nn ~= n
            error('sizes of x0 and A are incompatible')
        end
        b = b(:);
        if length(b) ~= m
            error('sizes of A and b are incompatible')
        end
    end

    % initialize multiplier variables
    if m > 0
        tauk = zeros(m,1);
    else
        tauk = [];  % return this
    end
    [tmp, g0] = f(xk);
    if size(g0,1) ~= n || size(g0,2) ~= 1
        error('gradient df(x) must be same size as x')
    end
    if all(g0 <= 0)
        mu0 = 1.0;
    else
        % if some components of g0 = grad f(x0) are positive
        % then average to generate a scale mu0; on average:  g0*x0 = mu0
        mu0 = sum(g0(g0 > 0) .* xk(g0 > 0)) / n;
    end
    lamk = mu0 ./ xk;

    % start output lists if requested
    if nargout > 3
        if m > 0
            iteratelist = [[xk; tauk; lamk]];
        else
            iteratelist = [[xk; lamk]];
        end
        nulist = [];
        mulist = [];
    end

    % loop: primal-dual interior point method uses Newton steps on barrier equations
    for k = 1:maxiters
        if any(xk <= 0.0)
            error('primal feasibility violated at iteration %d',k)
        end
        if any(lamk <= 0.0)
            error('dual feasibility violated at iteration %d',k)
        end
        [fk, gk, Hk] = f(xk);          % call user's function
        % check stopping condition
        if m > 0
            meritk = merit(xk,tauk,lamk,gk,A,b);
        else
            meritk = merit(xk,tauk,lamk,gk,[],[]);
        end
        if nargout > 3
            nulist = [nulist meritk];
        end
        if meritk < atol
            break
        end
        if k == 1
            merit0 = meritk;
        else
            if meritk / merit0 < rtol
                break
            end
        end
        % set up and solve Newton step equations
        mu = min(theta*meritk,meritk^2);
        if m > 0
            MM = [Hk,          -A',         -eye(n,n);
                  -A,          zeros(m,m),  zeros(m,n),
                  diag(lamk),  zeros(n,m),  diag(xk)];
            c = [- gk + A' * tauk + lamk;
                 A * xk - b;
                 mu - lamk .* xk];
        else
            MM = [Hk,         -eye(n,n);
                  diag(lamk), diag(xk)];
            c = [- gk + lamk;
                 mu - lamk .* xk];
        end
        %cond(MM)
        p = MM \ c;                          % Gaussian elimination: O(n^3)
        % apply ratio tests separately on x,lam
        % note: (dx,dtau,dlam) = (p(1:n),p(n+1:n+m),p(n+m+1:n+m+n))
        kappa = max(kappabar,1.0-meritk);
        [alphax, alphatau, alphalam] = ratiotest(xk,p(1:n),...
                                                 lamk,p(n+m+1:n+m+n),kappa);
        xk = xk + alphax * p(1:n);
        if m > 0
            tauk = tauk + alphatau * p(n+1:n+m);
        end
        lamk = lamk + alphalam * p(n+m+1:n+m+n);
        % append to output lists if desired
        if nargout > 3
            if m > 0
                iteratelist = [iteratelist, [xk; tauk; lamk]];
            else
                iteratelist = [iteratelist, [xk; lamk]];
            end
            % nulist already appended
            mulist = [mulist mu];
        end
    end
end

    function z = merit(x,tau,lam,dfx,A,b)
        if length(tau) > 0
            z = max([norm(dfx - A'*tau - lam),...
                     norm(b - A*x),...
                     norm(lam.*x)]);
        else
            z = max([norm(dfx - lam),...
                     norm(lam.*x)]);
        end
    end

    function [alphax, alphatau, alphalam] = ratiotest(x,dx,lam,dlam,kappa)
        alphax = 1.0;
        for j = 1:length(x)
            if dx(j) < 0
                alphax = min(alphax, - kappa * x(j) / dx(j));
            end
        end
        alphatau = 1.0;
        alphalam = 1.0;
        for j = 1:length(lam)
            if dlam(j) < 0
                alphalam = min(alphalam, - kappa * lam(j) / dlam(j));
            end
        end
    end
