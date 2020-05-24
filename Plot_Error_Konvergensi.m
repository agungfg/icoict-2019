clear all;
close all;
clc;

filename = 'TimeAndError.xlsx';
sheet = 4;

%Nx
range = 'B12:B18';
Nx = xlsread(filename,sheet,range);

%Mach
%range = 'B21:B26';
%L1_Mach = xlsread(filename,sheet,range);
range = 'F12:F18';
L2_Mach = xlsread(filename,sheet,range);

%Density
%range = 'C21:C26';
%L1_Density = xlsread(filename,sheet,range);
range = 'G12:G18';
L2_Density = xlsread(filename,sheet,range);


%plot(dx,L1_Mach,':','DisplayName','{\it L}^{1} Mach Number','color','black','linewidth',1.25);
plot(Nx,L2_Mach,'-','DisplayName','{\it L}^{2} Mach Number','color','black','linewidth',1.25);
hold on
%plot(dx,L1_Density,'-.','DisplayName','{\it L}^{1} Density','color','black','linewidth',1.25);
plot(Nx,L2_Density,'-.','DisplayName','{\it L}^{2} Density','color','black','linewidth',1.25);
axis([Nx(1) Nx(7) -7 -2])
xlabel('{\it Log (Nx)}')
ylabel('{\it Log (Error)}')
legend('show');
%legend('Location','southwest');