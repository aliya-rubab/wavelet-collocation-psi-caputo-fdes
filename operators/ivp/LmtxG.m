%==========================================================================
% LmtxG.m
% Author: Aliya Rubab
% Constructs the operational matrix associated with ?-Gegenbauer wavelets.
%
% This function evaluates the ?-Gegenbauer wavelet basis functions at the
% collocation points and assembles the corresponding operational matrix
% used in the wavelet collocation method for solving ?-Caputo fractional
% differential equations.
%
%Inputs:
%   k     - Resolution level
%   Len   - Length of the computational interval
%   M     - Number of wavelet basis functions
%   alp   - Gegenbauer parameter
%
% Note:
%   The input parameter alpha is retained for compatibility with the
%   overall implementation but is not used directly in this function.
%==========================================================================

function L=LmtxG(k, Len, alpha, M, alp)
%alpha=2;Len=1;k=2;M=5; alp=17.5;
N=2.^(k-1);
m=M*N;
f=@(x) x.^2;
for i=1:m
    x=Len*(2*i-1)./(2.*m);
    for n=1:N
         a=Len*(2.*n-2)./2.^k;
         b=Len*(2.*n)./2.^k;
        if x >= a && x < b
            if M==3
   FF(n,i)=(2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))+0*n;
   FF(n+N,i)=(2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1));
   FF(n+2*N,i)=(2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp));
          
           elseif M==4
   FF(n,i)=2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2))+0*n;
   FF(n+N,i)=2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1);
   FF(n+2*N,i)=2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp);  
   FF(n+3*N,i)=2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2);
            elseif M==5
   FF(n,i)=2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2))+0*n;
   FF(n+N,i)=2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1);
   FF(n+2*N,i)=2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp);  
   FF(n+3*N,i)=2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k.*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2);
   FF(n+4*N,i)=2.^(k./2).*(((96+24*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(4+2*alp))^(1/2))*(4*(2.^k.*f(x)-2*n+1)^4*alp+(22/3)*(2.^k.*f(x)-2*n+1)^4*alp^2-4*(2.^k.*f(x)-2*n+1)^2*alp+4*(2.^k.*f(x)-2*n+1)^4*alp^3-6*(2.^k.*f(x)-2*n+1)^2*alp^2+(2/3)*(2.^k.*f(x)-2*n+1)^4*alp^4-2*(2.^k.*f(x)-2*n+1)^2*alp^3+(1/2)*alp+(1/2)*alp^2);
        elseif M==6
   FF(n,i)=2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2))+0*n;
   FF(n+N,i)=2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1);
   FF(n+2*N,i)=2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp);  
   FF(n+3*N,i)=2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k.*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2);
   FF(n+4*N,i)=2.^(k./2).*(((96+24*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(4+2*alp))^(1/2))*(4*(2.^k.*f(x)-2*n+1)^4*alp+(22/3)*(2.^k.*f(x)-2*n+1)^4*alp^2-4*(2.^k.*f(x)-2*n+1)^2*alp+4*(2.^k.*f(x)-2*n+1)^4*alp^3-6*(2.^k.*f(x)-2*n+1)^2*alp^2+(2/3)*(2.^k.*f(x)-2*n+1)^4*alp^4-2*(2.^k.*f(x)-2*n+1)^2*alp^3+(1/2)*alp+(1/2)*alp^2);
   FF(n+5*N,i)=2.^(k./2).*(((600+120*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(5+2*alp))^(1/2))*((32/5)*(2.^k.*f(x)-2*n+1)^5*alp+(40/3)*(2.^k.*f(x)-2*n+1)^5*alp^2-8*(2.^k.*f(x)-2*n+1)^3*alp+(28/3)*(2.^k.*f(x)-2*n+1)^5*alp^3-(44/3)*(2.^k.*f(x)-2*n+1)^3*alp^2+(8/3)*(2.^k.*f(x)-2*n+1)^5*alp^4-8*(2.^k.*f(x)-2*n+1)^3*alp^3+2*alp*(2.^k.*f(x)-2*n+1)+3*(2.^k.*f(x)-2*n+1)*alp^2+(4/15)*(2.^k.*f(x)-2*n+1)^5*alp^5-(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^4+(2.^k.*f(x)-2*n+1)*alp^3);
          end
        elseif x > b
          FF(n,i)=0;  
        end
    end
end
L=FF;
