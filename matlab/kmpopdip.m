function kmpopdip(m,x0,rtol)
% KMPOPDIP  Solves the m-variable version of the Klee-Minty cube problem
% using POPDIP.  Compare
%   https://bueler.github.io/opt/assets/codes/kleeminty.m
% for a simplex method solution of this problem.  On this example, POPDIP
% seems to have exponentially-growing numbers of steps just like simplex.

    if nargin < 1,  m = 4;  end
    if nargin < 2,  x0 = ones(m,1);  end
    if nargin < 3,  rtol = 1.0e-12;  end

    % problem in "EZ" form:  Ax <= b
    A = eye(m,m);
    for k = 2:m
       A(k,1:k-1) = 2.^(k:-1:2);
    end
    b = 5.^(1:m)';
    % convert to standard form; now n = 2m
    A = [A eye(m,m)];
    x0 = [x0; b];

    % solve
    [xk,tauk,lamk,iterlist,nuklist,muklist] = popdip(...
        x0,@kleemintyfcn,A,b,rtol,1.0e-50,10000);

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
    format short g
    xk
    [fx,dfx,Hfx] = kleemintyfcn(xk);
    fx
end

    function [f,df,Hf] = kleemintyfcn(x)
        n = length(x);
        m = n / 2;
        c = -2.^(m-1:-1:0)';
        c = [c; zeros(m,1)];
        f = c' * x;
        df = c;
        Hf = zeros(n,n);
    end
