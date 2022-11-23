function obstacle(n,rtol)
% OBSTACLE  Uses POPDIP to solve a continuum optimizaton problem
% with positivity constraints:
%     min   f(u)
%     s.t.  u >= 0
% The quadratic function f(u) is given in OBSTACLEFCN.  This driver program
% calls POPDIP, prints information about the run, and generates two figures.

    if nargin < 1,  n = 20;  end
    if nargin < 2,  rtol = 1.0e-10;  end

    dx = 1/(n+1);
    x = dx:dx:1-dx;

    % initial iterate is strictly feasible
    u0 = 0.1 * ones(n,1);

    % solve using a faster-shrinking barrier: theta=0.001
    theta = 0.001;
    [uk,_,lamk,iterlist,nuklist,muklist] = popdip(u0,@obstaclefcn,[],[],rtol,1.0e-50,100,theta);

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
    fprintf('%d iterations with n=%d: ||u-uexact||_inf = %.3e\n',...
            size(iterlist,2),n,norm(uk-uex,'inf'))

    % plot iterates
    figure(1), clf
    K = size(iterlist,2);
    for j = 1:K-1
        plot(x,iterlist(1:n,j),'k'), hold on
    end
    plot(x,uk,'r'),  xlabel x,  grid on
    title('iterates in black, numerical solution in red')

    % plot exact solution and final iterate
    figure(2), clf
    subplot(3,1,1:2)
    xf = 0:.001:1;
    plot(x,uk,'r'),  hold on,  plot(xf,uexact(xf),'b')
    ylabel('solutions'),  grid on
    title('numerical solution in red, exact solution in blue')
    subplot(3,1,3)
    plot(x,abs(uk' - uexact(x)),'k')
    xlabel x,  ylabel('numerical error'),  grid on
end

    function uu = uexact(x)
        alf = 0.275562026630539;  % solves:  sin(2 pi alf)/(2 pi) + 0.7 (alf - 0.5)
        C = 100 * (cos(2*pi*alf) / (2*pi)^2 - (0.7/2) * alf * (alf - 1));
        uu = 100 * (-cos(2*pi*x) / (2*pi)^2 + (0.7/2) * x .* (x - 1)) + C;
        uu(x <= alf) = 0;
        uu(x >= 1-alf) = 0;
    end
