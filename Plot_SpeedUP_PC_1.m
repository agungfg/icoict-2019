clear all;
close all;
clc;
filename = 'TimeAndError.xlsx';
sheet = 1;
%Nx
range = 'A3:A9';
Nx = xlsread(filename,sheet,range);
%SpeedUP PC 1
range = 'H3:J9';
SP_1 = xlsread(filename,sheet,range);
plot(Nx,SP_1(:,1),'-','DisplayName','2 Threads','color','black','linewidth',1.25);
hold on
plot(Nx,SP_1(:,2),'--','DisplayName','4 Threads','color','black','linewidth',1.25);
plot(Nx,SP_1(:,3),'-.','DisplayName','8 Threads','color','black','linewidth',1.25);
axis([Nx(1) Nx(7) 0 8])
ylabel('SpeedUP')
xlabel('Number of grids (Nx)')
legend('show');
legend('Location','northwest');