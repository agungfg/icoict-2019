clear all;
close all;
clc;

%Initial Condition
dx      = 0.1;
%dx      = 0.00078125*2;
m=0.1/dx;
nFinal = 1401*m;
C       = 0.5; %Courant number
tFinal  = 20;
L       = 3;
N       = ( L/dx ) +1;

%Calculate Initial Condition
for i=1:N
    x(i)                = (i-1)*dx;
    P(i,1)              = 1-0.3146*x(i);
    T(i,1)              = 1-0.2314*x(i);
    V(i,1)              = (0.1+1.09*x(i))*T(i,1)^(1/2);
    A(i)                = 1 + 2.2*(x(i)-1.5)^2;
end;

t=0;
n=1;
while (t<tFinal)
%while (n<=nFinal)
    for i=2:N-1
        P_Pre(i)            =  (( -P(i,n)*( V(i+1,n) - V(i,n) ) ) / dx ) - ( P(i,n)*V(i,n)*(log(A(i+1))-log(A(i)))/dx) - (V(i,n)*(P(i+1,n)-P(i,n))/dx); 
        V_Pre(i)            =  (( -V(i,n)*( V(i+1,n) - V(i,n) ) ) / dx ) - ( (1/1.4) * ( ( ( T(i+1,n)- T(i,n)) / dx ) + (T(i,n)/P(i,n))*(P(i+1,n)-P(i,n))/dx)); 
        T_Pre(i)            =  (( -V(i,n)*( T(i+1,n) - T(i,n) ) ) / dx ) - ( (1.4-1)*T(i,n) * ( ( ( V(i+1,n)- V(i,n)) / dx ) + (V(i,n)*(log(A(i+1))-log(A(i)))/dx)));
        deltaT(i-1)         = C * (dx/((T(i,n)^(1/2))+V(i,n)));
    end
    minimunDeltaT = min(deltaT);
    for i=2:N-1
        P_FinalPre(i)       = P(i,n) + P_Pre(i)*minimunDeltaT;
        V_FinalPre(i)       = V(i,n) + V_Pre(i)*minimunDeltaT;
        T_FinalPre(i)       = T(i,n) + T_Pre(i)*minimunDeltaT;
    end
    P_FinalPre(1)           = 1;
    V_FinalPre(1)           = 2*V_FinalPre(2)-V_FinalPre(3);
    T_FinalPre(1)           = 1;
    for i=2:N-1
        P_Cor(i)            =  (( -P_FinalPre(i)*( V_FinalPre(i) - V_FinalPre(i-1) ) ) / dx ) - ( P_FinalPre(i)*V_FinalPre(i)*(log(A(i))-log(A(i-1)))/dx) - (V_FinalPre(i)*(P_FinalPre(i)-P_FinalPre(i-1))/dx);
        V_Cor(i)            =  (( -V_FinalPre(i)*( V_FinalPre(i) - V_FinalPre(i-1) ) ) / dx ) - ( (1/1.4) * ( ( ( T_FinalPre(i)- T_FinalPre(i-1)) / dx ) + (T_FinalPre(i)/P_FinalPre(i))*(P_FinalPre(i)-P_FinalPre(i-1))/dx)); 
        T_Cor(i)            =  (( -V_FinalPre(i)*( T_FinalPre(i) - T_FinalPre(i-1) ) ) / dx ) - ( (1.4-1)*T_FinalPre(i) * ( ( ( V_FinalPre(i)- V_FinalPre(i-1)) / dx ) + (V_FinalPre(i)*(log(A(i))-log(A(i-1)))/dx)));
    end
    for i=2:N-1
        P(i,n+1)            =  P(i,n) + 0.5*(P_Pre(i)+P_Cor(i))*minimunDeltaT;
        V(i,n+1)            =  V(i,n) + 0.5*(V_Pre(i)+V_Cor(i))*minimunDeltaT;
        T(i,n+1)            =  T(i,n) + 0.5*(T_Pre(i)+T_Cor(i))*minimunDeltaT;
    end
    
    %Left Boundary
    P(1,n+1)      = 1;
    V(1,n+1)      = 2*V(2,n+1)-V(3,n+1);
    T(1,n+1)      = 1;
    %Right Boundary
    P(N,n+1)      = 2*P(N-1,n+1)-P(N-2,n+1);
    V(N,n+1)      = 2*V(N-1,n+1)-V(N-2,n+1);
    T(N,n+1)      = 2*T(N-1,n+1)-T(N-2,n+1);
    
    %Plot Step
    t=t+minimunDeltaT;
    n=n+1;
    %plot(x,P(:,n),x,T(:,n),x,V(:,n))
    %plot(x,P(:,n),'DisplayName','Pressure','color','red','linewidth',1);
    %hold on
    %plot(x,T(:,n),'DisplayName','Temperature','color','green','linewidth',1);
    %plot(x,V(:,n),'DisplayName','Velocity','color','blue','linewidth',1);
    %legend('show');
    %title("Time at t = "+ num2str(t))
    %axis([0 3 0 2.5])
    %pause(0.00001);
    %hold off  
    
end

 %Plot Step
 %plot(x,P(:,n),x,T(:,n),x,V(:,n))
 plot(x,T(:,n),':','DisplayName','\itT','color','black','linewidth',1.25);
 hold on
 plot(x,V(:,n),'-.','DisplayName','\itV','color','black','linewidth',1.25);
 plot(x,P(:,n),'--','DisplayName','\rho','color','black','linewidth',1.25);
 legend('show');
 title("Time at t = "+ num2str(t))
 axis([0 3 0 2.5])
  
 %hold off  
 
