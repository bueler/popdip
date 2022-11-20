function linear(x0,rtol)
% LINEAR  Solves the small, 2-variable linear programming problem,
% in standard form, using the interior point method POPDIP:
%   min   [-1 -2 0 0 0] * x
%   s.t.  -2 x1 +   x2 + x3           = 2
%         -1 x1 + 2 x2 +    + x4      = 7
%            x1                  + x5 = 3
%         x >= 0
% This problem is solved by the simplex method in section 5.2 of
% Griva, Nash, Sofer (2009).

    % o.k. to be infeasible for linear eq. constraints:
    if nargin < 1,  x0 = [1 1 1 1 1];  end
    if nargin < 2,  rtol = 1.0e-14;  end

    A = [-2 1 1 0 0;
         -1 2 0 1 0;
          1 0 0 0 1];
    b = [2 7 3]';
    [xk,tauk,lamk,iterlist,nuklist,muklist] = popdip(x0,@linearfcn,A,b,rtol);

    %iterlist
    N = size(iterlist,2);
    x1 = iterlist(1,:);  x2 = iterlist(2,:);
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

    figure(1),  clf,  plot([0 0 1 3 3 0],[0 2 4 5 0 0],'r')
    hold on,  plot(x1,x2,'-ko'),  hold off
    axis([-1 4 -1 6]),  axis equal,  grid on
    xlabel('x_1','fontsize',20),  ylabel('x_2','fontsize',20)
end

    function [f,df,Hf] = linearfcn(x)
        c = [-1 -2 0 0 0]';
        f = c' * x;
        df = c;
        Hf = zeros(5,5);
    end
