function obstacle(n,rtol)
% OBSTACLE  Uses POPDIP to solve
%     min   f(u)
%     s.t.  u >= 0
% for f(u) given in OBSTACLEFCN.
% See OBSTACLEFCN and POPDIP.

    if nargin < 1,  n = 20;  end
    if nargin < 2,  rtol = 1.0e-12;  end

    dx = 1/(n+1);
    x = dx:dx:1-dx;

    % initial iterate is strictly feasible
    u0 = 0.1 * ones(n,1);

    % solve using a faster-shrinking barrier: theta=0.001
    theta = 0.001;
    [uk,_,lamk,iterlist,nuklist,muklist] = popdip(u0,@obstaclefcn,[],[],rtol,1.0e-50,1000,theta);

    % print run info
    N = size(iterlist,2);
    fprintf('        nu_k                 mu_k\n');
    for k = 1:N
        if k == 1
            fprintf('%3d: %20.15f\n',...
                    k-1,nuklist(k));
        else
            fprintf('%3d: %20.15f %20.15f\n',...
                    k-1,nuklist(k),muklist(k-1));
        end
    end

    % print summary and error norm
    uex = uexact(x)';
    fprintf('%d iterations with n=%d: ||u-uexact|| = %.3e\n',...
            size(iterlist,2),n,norm(uk-uex,'inf'))

    % plot iterates
    figure(1), clf
    subplot(1,2,1)
    K = size(iterlist,2);
    for j = 1:K-1
        plot(x,iterlist(1:n,j),'k'), hold on
    end
    plot(x,uk,'r')
    xlabel x, grid on

    % plot exact solution and final iterate
    subplot(1,2,2)
    plot(x,uk,'r')
    hold on
    xf = 0:.001:1;
    plot(xf,uexact(xf),'b')
    xlabel x, grid on
end

    function uu = uexact(x)
        a = 0.224437973369461;  % solve  sin(2 pi a) = 0.7 (2 pi a)
        C2 = (1/(2*pi)^2) * cos(pi*(1+2*a)) - (0.7/2) * (a + 0.5) * (a - 0.5);
        uu = -(1/(2*pi)^2) * cos(2*pi*x) + (0.7/2) * x .* (x - 1) + C2;
        uu(abs(x-0.5) >= a) = 0;
        uu = 100 * uu;
    end
