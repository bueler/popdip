function [f,df,Hf] = obstaclefcn(u)
% OBSTACLEFCN  Evaluate objective function, gradient, and Hessian for 1D
% obstacle problem.  The continuum objective is the elastic energy functional
%          /1
%   f[u] = |  (1/2) (u'(x))^2 - q(x) u(x) dx
%          /0
% Here q(x) = - 100 * (cos(2*pi*x) + 0.7), which is only positive in a small
% interval near x=0.5.  Uses midpoint rule.

    q = @(x) - 100 * (cos(2*pi*x) + 0.7);
    u = u(:);
    n = length(u);
    dx = 1.0 / (n+1);
    x = dx:dx:1.0-dx;   % length n

    % objective value
    f = 0;
    for i = 0:n         % loop over n+1 subintervals
        if i == 0
            dudx = (u(1) - 0)/dx;
            uav = 0.5 * (0 + u(1));
            xmid = 0.5 * (0 + x(1));
        elseif i == n
            dudx = (0 - u(n))/dx;
            uav = 0.5 * (u(n) + 0);
            xmid = 0.5 * (x(n) + 1);
        else
            dudx = (u(i+1) - u(i))/dx;
            uav = 0.5 * (u(i) + u(i+1));
            xmid = 0.5 * (x(i) + x(i+1));
        end
        f = f + (1/2) * dudx^2 - q(xmid) * uav;
    end
    f = f * dx;

    % gradient
    df = zeros(size(u));
    for i = 1:n         % loop over interior points
        if i == 1
            df(1) = (1/dx) * (-u(2) + 2*u(1));
        elseif i == n
            df(n) = (1/dx) * (2*u(n) - u(n-1));
        else
            df(i) = (1/dx) * (-u(i+1) + 2*u(i) - u(i-1));
        end
        df(i) = df(i) - (dx/2) * (q(x(i) - dx/2) + q(x(i) + dx/2));
    end

    % Hessian
    ee = ones(n,1);
    Hf = (1/dx) * full(spdiags([-ee, 2*ee, -ee],[-1,0,1],n,n));
end