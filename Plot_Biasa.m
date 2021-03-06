clear all;
close all;
clc;

%Initial Condition
dx      = 0.1;
C       = 0.5; %Courant number
Tfinal  = 28;
L       = 3;
N       = ( L/dx ) +1;

%Calculate Initial Condition
for i=1:N
    x(i)        = (i-1)*dx;
    P(i,1)      = 1-0.3146*x(i);
    T(i,1)      = 1-0.2314*x(i);
    V(i,1)      = (0.1+1.09*x(i))*T(i,1)^(1/2);
    A(i)        = 1 + 2.2*(x(i)-1.5)^2;
    P_Numerik(i,1)=P(i,1)*T(i,1);
end;

t=0;
n=1;
while (n<=1401)
    %First Predictor
    for i=2:N-1
        %Rho Scheme
        P_Pre(i,n+1)      =  (( -P(i,n)*( V(i+1,n) - V(i,n) ) ) / dx ) - ( P(i,n)*V(i,n)*(log(A(i+1))-log(A(i)))/dx) - (V(i,n)*(P(i+1,n)-P(i,n))/dx); 
        %V Scheme
        V_Pre(i,n+1)      =  (( -V(i,n)*( V(i+1,n) - V(i,n) ) ) / dx ) - ( (1/1.4) * ( ( ( T(i+1,n)- T(i,n)) / dx ) + (T(i,n)/P(i,n))*(P(i+1,n)-P(i,n))/dx)); 
        %T Scheme
        T_Pre(i,n+1)      =  (( -V(i,n)*( T(i+1,n) - T(i,n) ) ) / dx ) - ( (1.4-1)*T(i,n) * ( ( ( V(i+1,n)- V(i,n)) / dx ) + (V(i,n)*(log(A(i+1))-log(A(i)))/dx)));
        deltaT(i-1)        = C * (dx/((T(i,n)^(1/2))+V(i,n)));
    end
    
    %Left Boundary
    P_Pre(1,n+1)      = 1;
    V_Pre(1,n+1)      = 2*V_Pre(2,n+1)-V_Pre(3,n+1);
    T_Pre(1,n+1)      = 1;
    %Right Boundary
    P_Pre(N,n+1)      = 2*P_Pre(N-1,n+1)-P_Pre(N-2,n+1);
    V_Pre(N,n+1)      = 2*V_Pre(N-1,n+1)-V_Pre(N-2,n+1);
    T_Pre(N,n+1)      = 2*T_Pre(N-1,n+1)-T_Pre(N-2,n+1);
    
    %Pre Predictor
    minimunDeltaT = min(deltaT);
    for i=2:N-1
        P_FinalPre(i,n+1)      = P(i,n) + P_Pre(i,n+1)*minimunDeltaT;
        V_FinalPre(i,n+1)      = V(i,n) + V_Pre(i,n+1)*minimunDeltaT;
        T_FinalPre(i,n+1)      = T(i,n) + T_Pre(i,n+1)*minimunDeltaT;
    end
   
    %Left Boundary
    P_FinalPre(1,n+1)      = 1;
    V_FinalPre(1,n+1)      = 2*V_FinalPre(2,n+1)-V_FinalPre(3,n+1);
    T_FinalPre(1,n+1)      = 1;
    %Right Boundary
    P_FinalPre(N,n+1)      = 2*P_FinalPre(N-1,n+1)-P_FinalPre(N-2,n+1);
    V_FinalPre(N,n+1)      = 2*V_FinalPre(N-1,n+1)-V_FinalPre(N-2,n+1);
    T_FinalPre(N,n+1)      = 2*T_FinalPre(N-1,n+1)-T_FinalPre(N-2,n+1);
    
    %Corrector
    for i=2:N-1
        P_Cor(i,n+1)      =  (( -P_FinalPre(i,n+1)*( V_FinalPre(i,n+1) - V_FinalPre(i-1,n+1) ) ) / dx ) - ( P_FinalPre(i,n+1)*V_FinalPre(i,n+1)*(log(A(i))-log(A(i-1)))/dx) - (V_FinalPre(i,n+1)*(P_FinalPre(i,n+1)-P_FinalPre(i-1,n+1))/dx);
        V_Cor(i,n+1)      =  (( -V_FinalPre(i,n+1)*( V_FinalPre(i,n+1) - V_FinalPre(i-1,n+1) ) ) / dx ) - ( (1/1.4) * ( ( ( T_FinalPre(i,n+1)- T_FinalPre(i-1,n+1)) / dx ) + (T_FinalPre(i,n+1)/P_FinalPre(i,n+1))*(P_FinalPre(i,n+1)-P_FinalPre(i-1,n+1))/dx)); 
        T_Cor(i,n+1)      =  (( -V_FinalPre(i,n+1)*( T_FinalPre(i,n+1) - T_FinalPre(i-1,n+1) ) ) / dx ) - ( (1.4-1)*T_FinalPre(i,n+1) * ( ( ( V_FinalPre(i,n+1)- V_FinalPre(i-1,n+1)) / dx ) + (V_FinalPre(i,n+1)*(log(A(i))-log(A(i-1)))/dx)));
    end
    
    %Left Boundary
    P_Cor(1,n+1)      = 1;
    V_Cor(1,n+1)      = 2*V_Cor(2,n+1)-V_Cor(3,n+1);
    T_Cor(1,n+1)      = 1;
    %Right Boundary
    P_Cor(N,n+1)      = 2*P_Cor(N-1,n+1)-P_Cor(N-2,n+1);
    V_Cor(N,n+1)      = 2*V_Cor(N-1,n+1)-V_Cor(N-2,n+1);
    T_Cor(N,n+1)      = 2*T_Cor(N-1,n+1)-T_Cor(N-2,n+1);
    
    for i=2:N-1
        %Rho Scheme
        P(i,n+1)      =  P(i,n) + 0.5*(P_Pre(i,n+1)+P_Cor(i,n+1))*minimunDeltaT;
        %V Scheme
        V(i,n+1)      =  V(i,n) + 0.5*(V_Pre(i,n+1)+V_Cor(i,n+1))*minimunDeltaT;
        %T Scheme
        T(i,n+1)      =  T(i,n) + 0.5*(T_Pre(i,n+1)+T_Cor(i,n+1))*minimunDeltaT;
    end
    
    %Left Boundary
    P(1,n+1)      = 1;
    V(1,n+1)      = 2*V(2,n+1)-V(3,n+1);
    T(1,n+1)      = 1;
    %Right Boundary
    P(N,n+1)      = 2*P(N-1,n+1)-P(N-2,n+1);
    V(N,n+1)      = 2*V(N-1,n+1)-V(N-2,n+1);
    T(N,n+1)      = 2*T(N-1,n+1)-T(N-2,n+1);
    
    %Calculate Mach Number
    for i=1:N
        P_Numerik(i,n+1)=P(i,n+1)*T(i,n+1);
    end;
       
    %Plot Step
    t=t+minimunDeltaT;
    n=n+1;
    
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

    plot(x,P(:,n),'DisplayName','Pressure','color','red','linewidth',1);
    hold on
    plot(x,T(:,n),'DisplayName','Temperature','color','green','linewidth',1);
    plot(x,V(:,n),'DisplayName','Velocity','color','blue','linewidth',1);
    legend('show');
    title("Time at t = "+ num2str(t))
    axis([0 3 0 2.5])
    %pause(0.00001);
    %hold off   
