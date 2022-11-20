function obstacle(n,rtol)
% OBSTACLE  Uses POPDIP to solve
%     min   f(u)
%     s.t.  u >= 0
% for f(u) given in OBSTACLEFCN.

    if nargin < 1,  n = 20;         end
    if nargin < 2,  rtol = 1.0e-12;  end

    % solve
    u0 = 0.1 * ones(n,1);  % strictly feasible
    dx = 1/(n+1);
    x = dx:dx:1-dx;
    %u0 = uexact(dx:dx:1-dx) + 0.001;
    [uk,_,lamk,iterlist] = popdip(u0,@obstaclefcn,[],[],rtol);
    %format long, uklist'
    fprintf('%d iterations\n',size(iterlist,2))
    uex = uexact(x)';
    fprintf('||u-uexact|| = %.3e\n',norm(uk-uex,'inf'))

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
