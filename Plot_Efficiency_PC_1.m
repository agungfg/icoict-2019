clear all;
close all;
clc;

filename = 'TimeAndError.xlsx';
sheet = 1;

%dx
range = 'A3:A9';
Nx = xlsread(filename,sheet,range);

%Efficiency PC 1
range = 'L3:N9';
EP_1 = xlsread(filename,sheet,range);

plot(Nx,EP_1(:,1),'-','DisplayName','2 Threads','color','black','linewidth',1.25);
hold on
plot(Nx,EP_1(:,2),'--','DisplayName','4 Threads','color','black','linewidth',1.25);
plot(Nx,EP_1(:,3),'-.','DisplayName','8 Threads','color','black','linewidth',1.25);
axis([Nx(1) Nx(7) 20 100]);
ylabel('Efficiency (%)')
xlabel('Number of grids (Nx)')
legend('show');
legend('Location','southeast');