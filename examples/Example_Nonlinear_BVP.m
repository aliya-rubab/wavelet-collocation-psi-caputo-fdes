%==========================================================================
% Example1_Nonlinear_BVP.m
% Author: Aliya Rubab
%
% Example 1:
% Solves a nonlinear Psi-Caputo fractional boundary value problem using the
% psi-Gegenbauer wavelet collocation method.
%
% This example reproduces the numerical results reported in the paper.
%
% Dependencies:
%   Operators/BVP/LmtxG.m
%   Operators/BVP/PLFalphaG.m
%   Operators/BVP/KLG.m
%==========================================================================

function Example_Nonlinear_BVP
format short e
Len   = 1;
alpha = 1.9;
gama  = 0.9;
beta  = 1.9;
alp   = 17.5;
f = @(x) tan(x*pi./4);

%%%%% Parameter sets for the four cases
kSet  = [2 3 4 5];
MSet  = [3 4 5 6];
nnSet = [2 3 4 5];

%%%%% Run all four cases
for caseNo = 1:length(kSet)
    
    %%%%% Case parameters
    k  = kSet(caseNo);
    M  = MSet(caseNo);
    nn = nnSet(caseNo);
    
    N = 2.^(k-1);
    m = M*N;
    %%%% Initial solution
    sol  = zeros(m,nn);
    dsol = zeros(m,nn);

    %%% Nonlinear iteration
    for ii = 2:nn

        %%% Initial guess
        sol(:,1) = 0;
        sol(:,1) = 1;

        %%% Collocation points and problem terms
        for i = 1:m

            x = Len*(2*i-1)./(2*m);

            Yexact(i) = (f(x)-f(0)).^(2*alpha);
            xcol(i)   = x;

            a(i) = x;
            c(i) = exp(x);

            e(i) = x + exp(x).*sol(i,ii-1);
            j(i) = exp(x).*dsol(i,ii-1);

            A(i)=((f(x)-f(0))).*(((f(1)-f(0)).^(2*alpha))./(f(1)-f(0)));
            B(i)=(gamma(2)/gamma(2-gama)).*((f(x)-f(0)).^(1-gama)).*((f(1)-f(0)).^(2*alpha)./(f(1)-f(0)));

            F(i) = ...
                (gamma(2*alpha+1)/gamma(2*alpha+1-beta)).* ...
                (f(x)-f(0)).^(2*alpha-beta) ...
                + a(i).*(gamma(2*alpha+1)/ ...
                gamma(2*alpha+1-gama)).* ...
                (f(x)-f(0)).^(2*alpha-gama) ...
                + c(i).*(f(x)-f(0)).^(2*alpha).* ...
                (gamma(2*alpha+1)/ ...
                gamma(2*alpha+1-gama)).* ...
                (f(x)-f(0)).^(2*alpha-gama) ...
                + c(i).*sol(i,ii-1).*dsol(i,ii-1) ...
                - e(i).*(gamma(2)/gamma(2-gama)).* ...
                (f(x)-f(0)).^(1-gama).* ...
                ((f(1)-f(0)).^(2*alpha)./(f(1)-f(0))) ...
                - j(i).*(f(x)-f(0)).* ...
                ((f(1)-f(0)).^(2*alpha)./(f(1)-f(0)));

        end

        %%%% Operational matrices
        de = diag(e);
        dj = diag(j);

        L    = LmtxG(k,Len,alpha,M,alp);
        PLG  = PLFalphaG(k,alpha,Len,alp,M);
        PLG1 = PLFalphaG(k,beta,Len,alp,M);
        PLG2 = PLFalphaG(k,beta-gama,Len,alp,M);

        
        [KL KLL K1 K2]=KLG(Len, k, M, alpha, alp, beta, e, j, gama, sol, ii);

        %%%%% Solve for the nonlinear iteration
        Li = inv(L + PLG2*de + PLG1*dj - K1 - K2);

        C = F*Li;

        yn  = C*PLG1 + A - C*KL;
        dyn = C*(PLG2-KLL) + B;

        sol(:,ii)  = yn;
        dsol(:,ii) = dyn;

        %%%%%%%%% Absolute error
        for i = 1:m
            err(i) = abs(Yexact(i)-yn(i));
        end

    end

    %%%%%%%%% Plot numerical solution and absolute error
    figure(caseNo)

    subplot(2,1,1)
    plot(xcol,Yexact,'k',xcol,yn,'o','MarkerFaceColor','r')
    legend('Exact Solution','Numerical Solution')
    title(['u=',num2str(nn),', k=',num2str(k), ...
        ' & M=',num2str(M)])
    xlabel('x')
    ylabel('V(x)')

    subplot(2,1,2)
    plot(xcol,err,'b')
    xlabel('x')
    ylabel('Absolute Error')

end

end