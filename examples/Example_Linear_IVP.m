%==========================================================================
% Example_Linear_IVP.m
% Author: Aliya Rubab
%
% Example:
% Solves a linear Psi-Caputo fractional initial value problem using the
% psi-Gegenbauer wavelet collocation method.
%
% This example reproduces the numerical results for different values
% of the fractional order alpha.
%
% Dependencies:
%   Operators/BVP/LmtxG.m
%   Operators/BVP/PLFalphaG.m
%==========================================================================

function Example_Linear_IVP
format short e

Len = 1;
k   = 4;
M   = 5;
alp = 17.5;

f = @(x) x.^2;

%%%% Fractional-order parameter sets
alphaSet = [1.6 1.4 1.3 1.2 1.1];

%%%% Run all cases
for caseNo = 1:length(alphaSet)

    alpha = alphaSet(caseNo);

    %%%% Number of collocation points
    N = 2.^(k-1);
    m = M*N;

    %%% Preallocation
    Yexact = zeros(1,m);
    xcol   = zeros(1,m);
    b1     = zeros(1,m);
    V      = zeros(1,m);

    %%%% Collocation points and problem terms
    for i = 1:m

        x = Len*(2*i-1)./(2*m);

        Yexact(i) = ((f(x)-f(0)).^4) ...
                  - 2*((f(x)-f(0)).^3) - 7;

        xcol(i) = x;

        b1(i) = cos(f(x)).*sin(f(x));

        V(i) = (gamma(5)/gamma(5-alpha)) ...
             .* ((f(x)-f(0)).^(4-alpha)) ...
             - 2.*(gamma(4)/gamma(4-alpha)) ...
             .* ((f(x)-f(0)).^(3-alpha)) ...
             + b1(i).*( ...
             (gamma(5)/gamma(4)) ...
             .* ((f(x)-f(0)).^3) ...
             - 2*(gamma(4)/gamma(3)) ...
             .* ((f(x)-f(0)).^2) );

    end

    %%%% Operational matrices
    db1  = diag(b1);

    L    = LmtxG(k,Len,alpha,M,alp);
    PLG  = PLFalphaG(k,alpha,Len,alp,M);
    PLG1 = PLFalphaG(k,alpha-1,Len,alp,M);

    %%%% Numerical solution
    Li = inv(L + PLG1*db1);
    C  = V*Li;

    yn = C*PLG - 7;

    %%%% Absolute error
    error = abs(Yexact - yn);

    %%%% Plot exact and numerical solutions
    figure(caseNo)

    subplot(2,1,1)
    plot(xcol,Yexact,'m', ...
         xcol,yn,'.b','LineWidth',1.5)

    legend(['Exact sol at \alpha=',num2str(alpha)], ...
           ['Num sol at \alpha=',num2str(alpha)])

    xlabel('x')
    ylabel('V(x)')
    title(['Comparison of solutions at \alpha=',num2str(alpha)])

    %%%% Plot absolute error
    subplot(2,1,2)
    plot(xcol,error,'b')

    xlabel('x')
    ylabel('Absolute Error')
    title(['Absolute Error at \alpha=',num2str(alpha)])

end

end