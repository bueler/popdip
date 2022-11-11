function small(x0,tol)
% SMALL  Solves the small, 2-variable optimization problem
%     min   f(x) = (1/2) (x_1-1)^2 + (1/2) (x_2+1)^2
%     s.t.  x >= 0
% using POPDIP.  Prints digits; convergence looks quadratic.

if nargin < 1,  x0 = [2; 2];  end
if nargin < 2,  tol = 1.0e-14;  end

[xk,lamk,xklist,lamklist] = popdip(x0,@smallfcn,tol);

N = size(xklist,2);
for k = 1:N
    fprintf('%3d: %20.15f %20.15f\n',k,xklist(1,k),xklist(2,k));
end

figure(1),  clf
plot(xklist(1,:),xklist(2,:),'-ko')
maxx = max(xklist(1,:));
maxy = max(xklist(2,:));
axis([0 1.1*maxx 0 1.1*maxy])
grid on
xlabel('x_1','fontsize',20)
ylabel('x_2','fontsize',20)
end

    function [f,df,Hf] = smallfcn(x)
    % SMALLFCN  Quadratic function.  The unconstrained min is [1; -1].
    f = 0.5 * (x(1)-1)^2 + 0.5 * (x(2)+1)^2;
    df = [x(1)-1;
          x(2)+1];
    Hf = [1, 0;
          0, 1];
    end

