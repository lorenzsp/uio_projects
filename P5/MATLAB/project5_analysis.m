%% solutions with N = 100 & dt =1e-3
% solution using crank nicolson
file = 'cn.txt'
cn = importdata(file, ',');
% solution using implicit mthod
file = 'implicit.txt'
im = importdata(file, ',');
% analytical solution
file = 'solution.txt'
s = importdata(file, ',');

for ii = 1:10:1000
    clf
    hold on
    plot(0:0.01:1, (s(ii,:)),'g')
    plot(0:0.01:1, (cn(ii,:)),'k.')
    plot(0:0.01:1, (im(ii,:)),'r.')
    
    axis([0 1 0 1.2]);
    %plot(100:10:1000,A(1,ii));
    pause(0.01)
    
end

%% N = 10
file = 'data_10.txt'
A = importdata(file, ',');

figure();
hold on;
plot(A(:,1),A(:,2),'r--x','Markersize',8,'LineWidth',0.5);
plot(A(:,1),A(:,3),'m--x','Markersize',8,'LineWidth',0.5);
plot(A(:,1),A(:,4),'b--x','Markersize',8,'LineWidth',0.5);
ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error $\epsilon$ at $t=5\times10^{-2}$ and N=10','Interpreter','LateX','Fontsize',18);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',18);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 15)
grid on

figure();
hold on;
plot(A(:,1),A(:,5),'r--x','Markersize',8,'LineWidth',0.5);
plot(A(:,1),A(:,6),'m--x','Markersize',8,'LineWidth',0.5);
plot(A(:,1),A(:,7),'b--x','Markersize',8,'LineWidth',0.5);
ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error $\epsilon$ at $t=20\times10^{-2}$ and N=10','Interpreter','LateX','Fontsize',18);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',18);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 15)
grid on

%% N = 100
file = 'data_100.txt'
A = importdata(file, ',');


figure();
hold on;

plot(A(1:2,1),A(1:2,2),'r--x','Markersize',8,'LineWidth',0.5);
plot(A(:,1),A(:,3),'m--x','Markersize',8,'LineWidth',0.5);
plot(A(:,1),A(:,4),'b--x','Markersize',8,'LineWidth',0.5);
ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error $\epsilon$ at $t=5\times10^{-2}$ and N=100','Interpreter','LateX','Fontsize',18);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',18);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 15)
grid on


figure();
hold on;
plot(A(:,1),A(:,5),'r--x','Markersize',8,'LineWidth',0.5);
plot(A(:,1),A(:,6),'m--x','Markersize',8,'LineWidth',0.5);
plot(A(:,1),A(:,7),'b--x','Markersize',8,'LineWidth',0.5);
ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error $\epsilon$ at $t=20\times10^{-2}$ and N=100','Interpreter','LateX','Fontsize',18);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',18);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 15)
grid on

%% solution 10
file = 'solution_10.txt';
A = importdata(file, ',');

x = 0:0.1:1;

figure();
hold on;
plot(x, A(1,:),'o','Markersize',10)
plot(x, A(2,:),'x','Markersize',10)
plot(x, (A(3,:)),'d','Markersize',10);
plot(x, (A(4,:)),'-');

ll=legend(gca,'show','explicit','implicit','crank nicolson','analytical');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Solutions at $t=5\times10^{-2}$  $\Delta t=1\times10^{-5}$ N=10','Interpreter','LateX','Fontsize',18);
xlabel('x','interpreter','latex','fontsize',18);
ylabel('u(x,t)','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 15);
grid on

figure();
hold on;
plot(x, A(5,:),'o','Markersize',10)
plot(x, A(6,:),'x','Markersize',10)
plot(x, (A(7,:)),'d','Markersize',10);
plot(x, (A(8,:)),'-');

ll=legend(gca,'show','explicit','implicit','crank nicolson','analytical');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Solutions at $t=20\times10^{-2}$  $\Delta t=1\times10^{-5}$ N=10','Interpreter','LateX','Fontsize',18);
xlabel('x','interpreter','latex','fontsize',18);
ylabel('u(x,t)','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 15);
grid on


%% solution 100
file = 'solution_100.txt';
A = importdata(file, ',');

x = 0:0.01:1;

figure();
hold on;
plot(x, A(1,:),'o','Markersize',8)
plot(x, A(2,:),'x','Markersize',8)
plot(x, (A(3,:)),'d','Markersize',8);
plot(x, (A(4,:)),'-');

ll=legend(gca,'show','explicit','implicit','crank nicolson','analytical');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Solutions at $t=5\times10^{-2}$  $\Delta t=1\times10^{-5}$ N=100','Interpreter','LateX','Fontsize',18);
xlabel('x','interpreter','latex','fontsize',18);
ylabel('u(x,t)','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 15);
grid on

figure();
hold on;
plot(x, A(5,:),'o','Markersize',8)
plot(x, A(6,:),'x','Markersize',8)
plot(x, (A(7,:)),'d','Markersize',8);
plot(x, (A(8,:)),'-');


ll=legend(gca,'show','explicit','implicit','crank nicolson','analytical');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Solutions at $t=20\times10^{-2}$  $\Delta t=1\times10^{-5}$ N=100','Interpreter','LateX','Fontsize',18);
xlabel('x','interpreter','latex','fontsize',18);
ylabel('u(x,t)','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 15);
grid on

%% Crank nicolson matrix
%hold on;
figure();
file = 'CN_matrix.txt';
A = importdata(file, ',');
dt = log10([ 0.05,  0.005, 0.0005, 5e-05]);
dx = log10([ 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625, 0.00078125]);
surf(dt, dx,A,'Facecolor','interp')
zlim([-10 0]);
title('Relative error $\epsilon$ of Crank-Nicholson method at $t=0.8$','Interpreter','LateX','Fontsize',19);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',16);
ylabel('$\log_{10}(\Delta x)$','interpreter','latex','fontsize',16);
zlabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter','latex');
%% Implicit matrix

file = 'IM_matrix.txt';
A = importdata(file, ',');
dt = log10([ 0.05,  0.005, 0.0005, 5e-05]);
dx = log10([ 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625, 0.00078125]);
surf(dt, dx,A,'Facecolor','interp')
zlim([-10 0]);

title('Relative error $\epsilon$ of Implicit method','Interpreter','LateX','Fontsize',19);
hold on
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',16);
ylabel('$\log_{10}(\Delta x)$','interpreter','latex','fontsize',16);
zlabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter','latex');
%% Explicit matrix
% hold on;
figure();
file = 'EX_matrix.txt';
A = importdata(file, ',');
dt = log10([5e-06, 2e-05, 8e-05, 0.00032, 0.00128]);
dx = log10([ 0.1, 0.05, 0.025, 0.0125, 0.00625]);
surf(dt, dx,A,'Facecolor','interp');

zlim([-10 0]);
title('Relative error $\epsilon$ of Explicit method at $t=0.8$','Interpreter','LateX','Fontsize',19);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',16);
ylabel('$\log_{10}(\Delta x)$','interpreter','latex','fontsize',16);
zlabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter','latex');
%% 2d solution
file = 'sol_2d_N100_dt1e-3.txt';
B = importdata(file, ',');
figure();

xx = 0:0.01:1;
for ii = 0:1:max(size(B));
    
    q = 101*ii;
    qq = (ii+1)*101;
    surf(xx',xx',B((q+1):(qq),1:101));
    zlim([0 1])
    hold on
    pause(0.1)
    clf
end
%%

figure();
ii = 1;
q = 101*ii;
qq = (ii+1)*101;
surf(xx',xx',B((q+1):(qq),1:101),'FaceAlpha',0,'EdgeColor','interp','LineStyle',':');
hold on;
ii = 40;
q = 101*ii;
qq = (ii+1)*101;
surf(xx',xx',B((q+1):(qq),1:101));
colormap(spring)


zlim([0 01]);
title({'Numerical solutions of 2D diffusion equation'; 'at $t=0$  and $t = 0.03$ N = 100 and $\Delta t = 1\times 10^{-3}$'},'Interpreter','LateX','Fontsize',19);
xlabel('$x$','interpreter','latex','fontsize',16);
ylabel('$y$','interpreter','latex','fontsize',16);
zlabel('$u(x,y,t)$','interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter','latex');

%% 2d solution
file = 'relative_error_2d.txt';
B = importdata(file, ',');
%
%
dt = [-2, -3, -4, -5];
dx = log10(1./(10:10:80));
surf(dt,dx, B,'Facecolor','interp')
%
hold on;
title('Relative error $\epsilon$ of 2D diffusion equation with N = 50','Interpreter','LateX','Fontsize',19);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',16);
ylabel('$\log_{10}(\Delta x)$','interpreter','latex','fontsize',16);
zlabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',16)
grid on
xt = get(gca, 'XTick');
set(gca, 'FontSize', 15);
set(gca,'TickLabelInterpreter','latex');
