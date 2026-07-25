
%==========================================================================
% KLG.m
% Author: Aliya Rubab
%
% Computes the boundary-condition matrices required by the 
% psi-Gegenbauer wavelet collocation method for nonlinear
% Psi-Caputo fractional boundary value problems.
%
% Inputs:
%   Len    - Length of the computational interval
%   k      - Resolution level
%   M      - Number of wavelet basis functions
%   alpha  - Fractional order
%   alp    - Gegenbauer parameter
%   beta   - Problem parameter
%   e,j    - Auxiliary vectors
%   gama   - Fractional derivative order used in boundary conditions
%   sol    - Solution vector (if required)
%   ii     - Iteration index (if required)
% Notes:
%   beta, sol, and ii are retained for compatibility with the main program,
%   although they are not used directly in this function.
%
% Outputs:
%   KL, KLL, K1, K2 - Boundary condition matrices
%==========================================================================

function [KL KLL K1 K2]=KLG(Len, k, M, alpha, alp, beta, e, j, gama, sol, ii) 
%Len=1; k=2; M=5; alpha=2; beta=1; alp=17.5;
N=2.^(k-1);
m=M*N;
syms x 
 f = @(x) tan(x*pi./4);
   s = diff(f(x));
    for n = 1:N
         a = Len*(2.*n-2)./2.^k;
         b = Len*(2.*n)./2.^k;
        
            if M == 3
    L0=matlabFunction(((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
        V(n)=integral(L0,a,b);   
   
    L1=matlabFunction(((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
         
       V(n+N)=integral(L1,a,b);  
    
    L2=matlabFunction(((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
         
         
          V(n+2*N)=integral(L2,a,b);

              elseif M == 4
               
   L0=matlabFunction(((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
          
    V(n)=integral(L0,a,b);
    L1=matlabFunction(((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
         
   V(n+N)=integral(L1,a,b);
    L2=matlabFunction(((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
          
    V(n+2*N)=integral(L2,a,b);
   L3=matlabFunction((2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2)...
        .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
    
     V(n+3*N)=integral(L3,a,b);
            elseif M == 5
                
   L0=matlabFunction(((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
          
    V(n)=integral(L0,a,b);
    L1=matlabFunction(((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
         
   V(n+N)=integral(L1,a,b);
    L2=matlabFunction(((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
          
    V(n+2*N)=integral(L2,a,b);
    L3=matlabFunction((2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k.*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2)...
        .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
    
    V(n+3*N)=integral(L3,a,b);
   L4=matlabFunction((2.^(k./2).*(((96+24*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(4+2*alp))^(1/2))*(4*(2.^k.*f(x)-2*n+1)^4*alp+(22/3)*(2.^k.*f(x)-2*n+1)^4*alp^2-4*(2.^k.*f(x)-2*n+1)^2*alp+4*(2.^k.*f(x)-2*n+1)^4*alp^3-6*(2.^k.*f(x)-2*n+1)^2*alp^2+(2/3)*(2.^k.*f(x)-2*n+1)^4*alp^4-2*(2.^k.*f(x)-2*n+1)^2*alp^3+(1/2)*alp+(1/2)*alp^2)...
       .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
    V(n+4*N)=integral(L4,a,b);
    
            elseif M == 6
                
   L0=matlabFunction(((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
          
    V(n)=integral(L0,a,b);
    L1=matlabFunction(((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
         
   V(n+N)=integral(L1,a,b);
    L2=matlabFunction(((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
          
    V(n+2*N)=integral(L2,a,b);
    L3=matlabFunction((2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k.*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2)...
        .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
    
    V(n+3*N)=integral(L3,a,b);
   L4=matlabFunction((2.^(k./2).*(((96+24*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(4+2*alp))^(1/2))*(4*(2.^k.*f(x)-2*n+1)^4*alp+(22/3)*(2.^k.*f(x)-2*n+1)^4*alp^2-4*(2.^k.*f(x)-2*n+1)^2*alp+4*(2.^k.*f(x)-2*n+1)^4*alp^3-6*(2.^k.*f(x)-2*n+1)^2*alp^2+(2/3)*(2.^k.*f(x)-2*n+1)^4*alp^4-2*(2.^k.*f(x)-2*n+1)^2*alp^3+(1/2)*alp+(1/2)*alp^2)...
       .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
    V(n+4*N)=integral(L4,a,b);
          
    L5=matlabFunction((2.^(k./2).*(((600+120*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(5+2*alp))^(1/2))*((32/5)*(2.^k.*f(x)-2*n+1)^5*alp+(40/3)*(2.^k.*f(x)-2*n+1)^5*alp^2-8*(2.^k.*f(x)-2*n+1)^3*alp+(28/3)*(2.^k.*f(x)-2*n+1)^5*alp^3-(44/3)*(2.^k.*f(x)-2*n+1)^3*alp^2+(8/3)*(2.^k.*f(x)-2*n+1)^5*alp^4-8*(2.^k.*f(x)-2*n+1)^3*alp^3+2*alp*(2.^k.*f(x)-2*n+1)+3*(2.^k.*f(x)-2*n+1)*alp^2+(4/15)*(2.^k.*f(x)-2*n+1)^5*alp^5-(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^4+(2.^k.*f(x)-2*n+1)*alp^3)...
       .*((f(1)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha)));
     V(n+5*N)=integral(L5,a,b);   
            end
    end
    % Construct boundary-condition vectors
for i = 1:m
   x = Len*(2.*i-1)./(2*m);


    B(1,i)=(f(x)-f(0))./(f(1)-f(0));
    BB(1,i)=(gamma(2)/gamma(2-gama)).*(((f(x)-f(0)).^(1-gama))./(f(1)-f(0)));
    B1(1,i)=j(i).*((f(x)-f(0)).^(1))./(f(1)-f(0));
    B2(1,i)=e(i).*(gamma(2)/gamma(2-gama)).*((f(x)-f(0)).^(1-gama))./(f(1)-f(0));
    
end
KL=V'*B;
KLL=V'*BB;
K1=V'*B1;
K2=V'*B2;
end
