function small(x0,rtol)
% SMALL  Solves the small, 2-variable quadratic optimization problem
%     min   f(x) = (1/2) (x_1-1)^2 + (1/2) (x_2+1)^2
%     s.t.  x >= 0
% using POPDIP.  Exact solution x* = [1; 0].  Prints digits.  Convergence
% looks quadratic.

    if nargin < 1,  x0 = [2; 2];  end
    if nargin < 2,  rtol = 1.0e-14;  end

    [xk,tauk,lamk,iterlist,nuklist,muklist] = popdip(x0,@smallfcn,[],[],rtol);

    N = size(iterlist,2);
    x1 = iterlist(1,:);  x2 = iterlist(2,:);
    lam1 = iterlist(3,:);  lam2 = iterlist(4,:);

    % print primal iterates, merif function values, and barrier parameters
    fprintf('        x_1                  x_2                  nu_k                 mu_k\n');
    for k = 1:N
        if k == 1
            fprintf('%3d: %20.15f %20.15f %20.15f\n',...
                    k-1,x1(k),x2(k),nuklist(k));
        else
            fprintf('%3d: %20.15f %20.15f %20.15f %20.15f\n',...
                    k-1,x1(k),x2(k),nuklist(k),muklist(k-1));
        end
    end

    % primal iterates in figure 1, dual in figure 2
    figure(1),  clf,  plot(x1,x2,'-ko')
    axis([0 1.1*max(x1) 0 1.1*max(x2)]),  grid on
    xlabel('x_1','fontsize',20),  ylabel('x_2','fontsize',20)
    figure(2),  clf,  plot(lam1,lam2,'-ko')
    axis([0 1.1*max(lam1) 0 1.1*max(lam2)]),  grid on
    xlabel('\lambda_1','fontsize',20),  ylabel('\lambda_2','fontsize',20)
end

    function [f,df,Hf] = smallfcn(x)
    % SMALLFCN  Quadratic function.  The unconstrained min is [1; -1].
        f = 0.5 * (x(1)-1)^2 + 0.5 * (x(2)+1)^2;
        df = [x(1)-1;
              x(2)+1];
        Hf = [1, 0;
              0, 1];
    end
