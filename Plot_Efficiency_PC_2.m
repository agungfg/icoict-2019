clear all;
close all;
clc;

filename = 'TimeAndError.xlsx';
sheet = 1;

%Nx
range = 'A3:A9';
Nx = xlsread(filename,sheet,range);

%Efficiency PC 2
range = 'L12:N18';
EP_2 = xlsread(filename,sheet,range);

plot(Nx,EP_2(:,1),'-','DisplayName','2 Threads','color','black','linewidth',1.25);
hold on
plot(Nx,EP_2(:,2),'--','DisplayName','4 Threads','color','black','linewidth',1.25);
plot(Nx,EP_2(:,3),'-.','DisplayName','8 Threads','color','black','linewidth',1.25);
axis([Nx(1) Nx(7) 20 100])
ylabel('Efficiency (%)')
xlabel('Number of grids (Nx)')
legend('show');
legend('Location','southeast');