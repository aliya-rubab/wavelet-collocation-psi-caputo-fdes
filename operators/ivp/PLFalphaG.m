%==========================================================================
% PLFalphaG.m
% Author: Aliya Rubab.
%
% Computes the psi-Caputo fractional integration operational matrix for
% psi-Gegenbauer wavelets.
%
% This function evaluates the psi-Caputo fractional integral of the
% psi-Gegenbauer wavelet basis functions using numerical quadrature and
% constructs the corresponding operational matrix used in the wavelet
% collocation method.
%
% Inputs:
%   k      - Resolution level
%   alpha  - Fractional order
%   Len    - Length of the computational interval
%   alp    - Gegenbauer parameter
%   M      - Number of wavelet basis functions
%
% Output:
%   PLG    - Operational matrix representing the psi-Caputo fractional
%            integral of the psi-Gegenbauer wavelet basis
%==========================================================================
%==========================================================================

function PLG=PLFalphaG(k,alpha,Len,alp,M)
%alpha=2;Len=1;k=2;M=5; alp=17.5;
N=2.^(k-1);
m=M*N;
 syms x 
  f=@(x) x.^2;
   s= diff(f(x));
for i=1:m
    t=Len*(2*i-1)./(2.*m);
    for n=1:N
         a=Len*(2.*n-2)./2.^k;
         b=Len*(2.*n)./2.^k;
        if t >= a && t < b
            if M==3
    L0=matlabFunction((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n)=integral(L0,a,t);
    L1=matlabFunction((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
         
    V(i,n+N)=integral(L1,a,t);
    L2=matlabFunction((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n+2.*N)=integral(L2,a,t);
            elseif M==4
               
   L0=matlabFunction((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n)=integral(L0,a,t);
    L1=matlabFunction((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
         
    V(i,n+N)=integral(L1,a,t);
    L2=matlabFunction((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n+2.*N)=integral(L2,a,t);
   L3=matlabFunction(2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2)...
        .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
    
     V(i,n+3.*N)=integral(L3,a,t);
            elseif M==5
                
   L0=matlabFunction((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n)=integral(L0,a,t);
    L1=matlabFunction((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
         
    V(i,n+N)=integral(L1,a,t);
    L2=matlabFunction((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n+2.*N)=integral(L2,a,t);
    L3=matlabFunction(2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k.*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2)...
        .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
    
    V(i,n+3.*N)=integral(L3,a,t);
   L4=matlabFunction(2.^(k./2).*(((96+24*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(4+2*alp))^(1/2))*(4*(2.^k.*f(x)-2*n+1)^4*alp+(22/3)*(2.^k.*f(x)-2*n+1)^4*alp^2-4*(2.^k.*f(x)-2*n+1)^2*alp+4*(2.^k.*f(x)-2*n+1)^4*alp^3-6*(2.^k.*f(x)-2*n+1)^2*alp^2+(2/3)*(2.^k.*f(x)-2*n+1)^4*alp^4-2*(2.^k.*f(x)-2*n+1)^2*alp^3+(1/2)*alp+(1/2)*alp^2)...
       .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
    V(i,n+4.*N)=integral(L4,a,t);
    
            elseif M==6
                
   L0=matlabFunction((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n)=integral(L0,a,t);
    L1=matlabFunction((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
         
    V(i,n+N)=integral(L1,a,t);
    L2=matlabFunction((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n+2.*N)=integral(L2,a,t);
    L3=matlabFunction(2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k.*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2)...
        .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
    
    V(i,n+3.*N)=integral(L3,a,t);
   L4=matlabFunction(2.^(k./2).*(((96+24*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(4+2*alp))^(1/2))*(4*(2.^k.*f(x)-2*n+1)^4*alp+(22/3)*(2.^k.*f(x)-2*n+1)^4*alp^2-4*(2.^k.*f(x)-2*n+1)^2*alp+4*(2.^k.*f(x)-2*n+1)^4*alp^3-6*(2.^k.*f(x)-2*n+1)^2*alp^2+(2/3)*(2.^k.*f(x)-2*n+1)^4*alp^4-2*(2.^k.*f(x)-2*n+1)^2*alp^3+(1/2)*alp+(1/2)*alp^2)...
       .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
    V(i,n+4.*N)=integral(L4,a,t);
          
    L5=matlabFunction(2.^(k./2).*(((600+120*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(5+2*alp))^(1/2))*((32/5)*(2.^k.*f(x)-2*n+1)^5*alp+(40/3)*(2.^k.*f(x)-2*n+1)^5*alp^2-8*(2.^k.*f(x)-2*n+1)^3*alp+(28/3)*(2.^k.*f(x)-2*n+1)^5*alp^3-(44/3)*(2.^k.*f(x)-2*n+1)^3*alp^2+(8/3)*(2.^k.*f(x)-2*n+1)^5*alp^4-8*(2.^k.*f(x)-2*n+1)^3*alp^3+2*alp*(2.^k.*f(x)-2*n+1)+3*(2.^k.*f(x)-2*n+1)*alp^2+(4/15)*(2.^k.*f(x)-2*n+1)^5*alp^5-(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^4+(2.^k.*f(x)-2*n+1)*alp^3)...
       .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
     V(i,n+5.*N)=integral(L5,a,t);   
            end
                   
   elseif t>b  
        if M==3

          L0=matlabFunction((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          V(i,n)=integral(L0,a,b);
          L1=matlabFunction((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          V(i,n+N)=integral(L1,a,b);
          L2=matlabFunction((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
             .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          V(i,n+2.*N)=integral(L2,a,b);
          
        elseif M==4
           
    L0=matlabFunction((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n)=integral(L0,a,b);
    L1=matlabFunction((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
         
    V(i,n+N)=integral(L1,a,b);
    L2=matlabFunction((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n+2.*N)=integral(L2,a,b);
   L3=matlabFunction(2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2)...
        .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
    
     V(i,n+3.*N)=integral(L3,a,b);
        
        elseif M==5
   L0=matlabFunction((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n)=integral(L0,a,b);
    L1=matlabFunction((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
         
    V(i,n+N)=integral(L1,a,b);
    L2=matlabFunction((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n+2.*N)=integral(L2,a,b);
    L3=matlabFunction(2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2)...
        .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
    
     V(i,n+3.*N)=integral(L3,a,b);
    L4=matlabFunction(2.^(k./2).*(((96+24*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(4+2*alp))^(1/2))*(4*(2.^k.*f(x)-2*n+1)^4*alp+(22/3)*(2.^k.*f(x)-2*n+1)^4*alp^2-4*(2.^k.*f(x)-2*n+1)^2*alp+4*(2.^k.*f(x)-2*n+1)^4*alp^3-6*(2.^k.*f(x)-2*n+1)^2*alp^2+(2/3)*(2.^k.*f(x)-2*n+1)^4*alp^4-2*(2.^k.*f(x)-2*n+1)^2*alp^3+(1/2)*alp+(1/2)*alp^2)...
       .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
     V(i,n+4.*N)=integral(L4,a,b);
         
        elseif M==6
                
    L0=matlabFunction((2.^(k./2).*((alp*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2*alp))^(1/2)))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n)=integral(L0,a,b);
    L1=matlabFunction((2.^(k./2).*(((alp+1)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(1+2*alp))^(1/2))*2*alp*(2.^k.*f(x)-2*n+1))...
             .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
         
    V(i,n+N)=integral(L1,a,b);
    L2=matlabFunction((2.^(k./2).*(((4+2*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(2+2*alp))^(1/2))*(2*alp*(1+alp).*(2.^k.*f(x)-2*n+1).^2-alp))...
              .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
          
    V(i,n+2.*N)=integral(L2,a,b);
    L3=matlabFunction(2.^(k./2).*(((18+6*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(3+2*alp))^(1/2))*((8/3)*(2.^k*f(x)-2*n+1)^3*alp+4*(2.^k.*f(x)-2*n+1)^3*alp^2-2*alp*(2.^k.*f(x)-2*n+1)+(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^3-2*(2.^k.*f(x)-2*n+1)*alp^2)...
        .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
    
     V(i,n+3.*N)=integral(L3,a,b);
    L4=matlabFunction(2.^(k./2).*(((96+24*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(4+2*alp))^(1/2))*(4*(2.^k.*f(x)-2*n+1)^4*alp+(22/3)*(2.^k.*f(x)-2*n+1)^4*alp^2-4*(2.^k.*f(x)-2*n+1)^2*alp+4*(2.^k.*f(x)-2*n+1)^4*alp^3-6*(2.^k.*f(x)-2*n+1)^2*alp^2+(2/3)*(2.^k.*f(x)-2*n+1)^4*alp^4-2*(2.^k.*f(x)-2*n+1)^2*alp^3+(1/2)*alp+(1/2)*alp^2)...
       .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
     V(i,n+4.*N)=integral(L4,a,b);
     
    L5=matlabFunction(2.^(k./2).*(((600+120*alp)*gamma(alp)^2/pi/(2^(1-2*alp))/gamma(5+2*alp))^(1/2))*((32/5)*(2.^k.*f(x)-2*n+1)^5*alp+(40/3)*(2.^k.*f(x)-2*n+1)^5*alp^2-8*(2.^k.*f(x)-2*n+1)^3*alp+(28/3)*(2.^k.*f(x)-2*n+1)^5*alp^3-(44/3)*(2.^k.*f(x)-2*n+1)^3*alp^2+(8/3)*(2.^k.*f(x)-2*n+1)^5*alp^4-8*(2.^k.*f(x)-2*n+1)^3*alp^3+2*alp*(2.^k.*f(x)-2*n+1)+3*(2.^k.*f(x)-2*n+1)*alp^2+(4/15)*(2.^k.*f(x)-2*n+1)^5*alp^5-(4/3)*(2.^k.*f(x)-2*n+1)^3*alp^4+(2.^k.*f(x)-2*n+1)*alp^3)...
       .*((f(t)-f(x)).^(alpha-1)).*subs(s,x,x)./gamma(alpha));
     V(i,n+5.*N)=integral(L5,a,b);
            end
        end
    end
end
PLG=V';
