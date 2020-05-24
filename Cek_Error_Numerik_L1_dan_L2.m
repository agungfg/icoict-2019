clear all;
close all;
clc;

%Initial Condition
N       = 50;
C       = 0.5; %Courant number
Tfinal  = 28.9516;
L       = 3;
dx      = (L/(N-1));

%Calculate Initial Condition
%parpool('local',8)
for i=1:N
    x(i)        = (i-1)*dx;
    P_Old(i)      = 1-0.3146*x(i);
    T_Old(i)      = 1-0.2314*x(i);
    V_Old(i)      = (0.1+1.09*x(i))*T_Old(i)^(1/2);
    A(i)        = 1 + 2.2*(x(i)-1.5)^2;
end;

t=0;
n=1;
while (t<Tfinal)
    %First Predictor
    for i=2:N-1
        %Rho Scheme
        P_Pre(i)      =  (( -P_Old(i)*( V_Old(i+1) - V_Old(i) ) ) / dx ) - ( P_Old(i)*V_Old(i)*(log(A(i+1))-log(A(i)))/dx) - (V_Old(i)*(P_Old(i+1)-P_Old(i))/dx); 
        %V Scheme
        V_Pre(i)      =  (( -V_Old(i)*( V_Old(i+1) - V_Old(i) ) ) / dx ) - ( (1/1.4) * ( ( ( T_Old(i+1)- T_Old(i)) / dx ) + (T_Old(i)/P_Old(i))*(P_Old(i+1)-P_Old(i))/dx)); 
        %T Scheme
        T_Pre(i)      =  (( -V_Old(i)*( T_Old(i+1) - T_Old(i) ) ) / dx ) - ( (1.4-1)*T_Old(i) * ( ( ( V_Old(i+1)- V_Old(i)) / dx ) + (V_Old(i)*(log(A(i+1))-log(A(i)))/dx)));
        deltaT(i-1)        = C * (dx/((T_Old(i)^(1/2))+V_Old(i)));
    end
    
    %Left Boundary
    P_Pre(1)      = 1;
    V_Pre(1)      = 2*V_Pre(2)-V_Pre(3);
    T_Pre(1)      = 1;
    %Right Boundary
    P_Pre(N)      = 2*P_Pre(N-1)-P_Pre(N-2);
    V_Pre(N)      = 2*V_Pre(N-1)-V_Pre(N-2);
    T_Pre(N)      = 2*T_Pre(N-1)-T_Pre(N-2);
    
    %Pre Predictor
    minimunDeltaT = min(deltaT);
    for i=2:N-1
        P_FinalPre(i)      = P_Old(i) + P_Pre(i)*minimunDeltaT;
        V_FinalPre(i)      = V_Old(i) + V_Pre(i)*minimunDeltaT;
        T_FinalPre(i)      = T_Old(i) + T_Pre(i)*minimunDeltaT;
    end
   
    %Left Boundary
    P_FinalPre(1)      = 1;
    V_FinalPre(1)      = 2*V_FinalPre(2)-V_FinalPre(3);
    T_FinalPre(1)      = 1;
    %Right Boundary
    P_FinalPre(N)      = 2*P_FinalPre(N-1)-P_FinalPre(N-2);
    V_FinalPre(N)      = 2*V_FinalPre(N-1)-V_FinalPre(N-2);
    T_FinalPre(N)      = 2*T_FinalPre(N-1)-T_FinalPre(N-2);
    
    %Corrector
    for i=2:N-1
        P_Cor(i)      =  (( -P_FinalPre(i)*( V_FinalPre(i) - V_FinalPre(i-1) ) ) / dx ) - ( P_FinalPre(i)*V_FinalPre(i)*(log(A(i))-log(A(i-1)))/dx) - (V_FinalPre(i)*(P_FinalPre(i)-P_FinalPre(i-1))/dx);
        V_Cor(i)      =  (( -V_FinalPre(i)*( V_FinalPre(i) - V_FinalPre(i-1) ) ) / dx ) - ( (1/1.4) * ( ( ( T_FinalPre(i)- T_FinalPre(i-1)) / dx ) + (T_FinalPre(i)/P_FinalPre(i))*(P_FinalPre(i)-P_FinalPre(i-1))/dx)); 
        T_Cor(i)      =  (( -V_FinalPre(i)*( T_FinalPre(i) - T_FinalPre(i-1) ) ) / dx ) - ( (1.4-1)*T_FinalPre(i) * ( ( ( V_FinalPre(i)- V_FinalPre(i-1)) / dx ) + (V_FinalPre(i)*(log(A(i))-log(A(i-1)))/dx)));
    end
    
    %Left Boundary
    P_Cor(1)      = 1;
    V_Cor(1)      = 2*V_Cor(2)-V_Cor(3);
    T_Cor(1)      = 1;
    %Right Boundary
    P_Cor(N)      = 2*P_Cor(N-1)-P_Cor(N-2);
    V_Cor(N)      = 2*V_Cor(N-1)-V_Cor(N-2);
    T_Cor(N)      = 2*T_Cor(N-1)-T_Cor(N-2);
    
    for i=2:N-1
        %Rho Scheme
        P_New(i)      =  P_Old(i) + 0.5*(P_Pre(i)+P_Cor(i))*minimunDeltaT;
        %V Scheme
        V_New(i)      =  V_Old(i) + 0.5*(V_Pre(i)+V_Cor(i))*minimunDeltaT;
        %T Scheme
        T_New(i)      =  T_Old(i) + 0.5*(T_Pre(i)+T_Cor(i))*minimunDeltaT;
    end
    
    %Left Boundary
    P_New(1)      = 1;
    V_New(1)      = 2*V_New(2)-V_New(3);
    T_New(1)      = 1;
    %Right Boundary
    P_New(N)      = 2*P_New(N-1)-P_New(N-2);
    V_New(N)      = 2*V_New(N-1)-V_New(N-2);
    T_New(N)      = 2*T_New(N-1)-T_New(N-2);
    
    %Plot Step
    t=t+minimunDeltaT;
    n=n+1;
    
    P_Old = P_New;
    T_Old = T_New;
    V_Old = V_New;
end
%hitung error
L1_Mach = 0;
L2_Mach = 0;
L1_Density= 0;
L2_Density = 0;
for i=1:N
    
    M_Numerik(i)=V_New(i)/sqrt(T_New(i));
    P_Numerik(i)=P_New(i)*T_New(i);
    
    syms m;
    eqn = 1+2.2*(x(i)-1.5).^2 ==((5+m.^2).^3)/(216*m);
    mach = solve(eqn,m);
    mach = double(mach);
    if (x(i) <= 1.5)
        M_Exact(i)=mach(1);
    else
        M_Exact(i)=mach(2);
    end;   
    P_Exact(i)=(1+((1.4-1)*M_Exact(i)^2)/2)^((-1)/(1.4-1));
    M_difference(i)=100*abs((M_Exact(i) - M_Numerik(i))/M_Exact(i));
    P_difference(i)=100*(abs(P_Exact(i) - P_New(i))/P_Exact(i));
    
    L1_Mach = L1_Mach + abs(M_Exact(i) - M_Numerik(i))*dx;
    L2_Mach = L2_Mach + ((abs(M_Exact(i) - M_Numerik(i)).^2)*dx);
    
    L1_Density = L1_Density + abs(P_Exact(i) - P_New(i))*dx;
    L2_Density = L2_Density + ((abs(P_Exact(i) - P_New(i)).^2)*dx);
  
end;

Convergences = [L1_Mach, L1_Density, L2_Mach.^0.5, L2_Density.^0.5];

fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

%title('Exact and Numerical Result')
xlabel('Grid \it x')

yyaxis left

plot(x,P_Exact,'*','color','black','linewidth',1.25);
hold on
plot(x,P_New,'color','black','linewidth',1.25);
ylabel('Density \rho / \rho_{0}')
axis([0 3 0 1.05]);

yyaxis right
ylabel('Mach Number')
plot(x,M_Exact,'*','color','black','linewidth',1.25);
hold on
plot(x,M_Numerik,'color','black','linewidth',1.25);
axis([0 3 0 3.6]);

legend('Exact','Numerical');
legend('Location','north')
