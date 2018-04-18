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
<<<<<<< HEAD
%%
% analytical solution
file = 'test.txt'
s = importdata(file, ',');

for ii = 1:1:1000
    clf
    hold on
    plot(0:0.01:1, s(ii,:),'g')
    %plot(0:0.01:1, (cn(ii,:)),'k.')
    %plot(0:0.01:1, (im(ii,:)),'r.')
    
    axis([0 1 0 1.2]);
    %plot(100:10:1000,A(1,ii));
    pause(0.01)
    
end

=======
>>>>>>> 1d184eeb16653cce57841c04082d42d2318dccf1

%% speranza 10
file = 'espiazione_ritorno10.txt'
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

%% speranza 100
file = 'espiazione_ritorno100.txt'
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

%% espiazione 10
file = 'espiazione10.txt';
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


%% espiazione 100
file = 'espiazione100.txt';
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

%%






















file = 'cnmerda.txt';
A = importdata(file, ',');

for ii =1:1:1000
    hold on
plot(A(1,:),A(ii,:));
pause(0.001)
clf
end


file = 'mado640.txt';
A = importdata(file, ',');

figure();
 hold on;
 
plot(A(:,1),A(:,2),'rx');
plot(A(:,1),A(:,3),'bx');
plot(A(:,1),A(:,4),'gx');
%plot(A(1,:),A(5,:),'-');
% error



figure();
hold on;
plot(dx:dx:1,abs(A(1,2:length())-A(4,2:641)./A(4,2:641)),'x');
plot(dx:dx:1,abs(A(2,2:641)-A(4,2:641)./A(4,2:641)),'x');
plot(dx:dx:1,abs(A(3,2:641)-A(4,2:641)./A(4,2:641)),'x');



%x = -1:-1:-7;

%plot(x, (A(4,:)),'-');

ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error $\Delta t=1e-4$ at $t=0.99$','Interpreter','LateX','Fontsize',18);
xlabel('x','interpreter','latex','fontsize',15);
ylabel('u(x,t)','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter','latex');
grid on



























%% per giovanni con affetto
%
%%
file = 'riprova.txt'
A = importdata(file, ',');

dt = A(:,1);

hold on;
plot(dt,A(:,2),'rx')
plot(dt, (A(:,3)),'bx');
plot(dt, (A(:,4)),'cx');

ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error as a function of $\Delta t$ at $t=5\times10^{-2}$ N=10','Interpreter','LateX','Fontsize',18);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',15);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter','latex');
grid on


%% dx = 1/10
% in x = 0.8
%% t = 50e-3
file = 'ex_im_cn_N10_50&200e-3.txt';
A = importdata(file, ',');

dt = log10(5e-3.*[1, 1./5, 1./25, 1./125, 1./625, 1./3125, 1./15625]);
figure();
hold on;
plot(dt,A(:,2),'r-')
plot(dt, (A(:,3)),'b-');
plot(dt, (A(:,4)),'c-');

ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error as a function of $\Delta t$ at $t=5\times10^{-2}$ N=10','Interpreter','LateX','Fontsize',18);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',15);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter','latex');
grid on


% file = 'ex_im_cn_an_t50e-3_dt5e-3_N10.txt';
% A = importdata(file, ',');
% figure();
% 
% x = 0:0.1:1;
% figure();
% hold on;
% plot(x, A(1,:),'-')
% plot(x, A(2,:),'-')
% plot(x, (A(3,:)),'-');
% plot(x, (A(4,:)),'-');
% 


%% t=200e-3
% sbagliato lexplicit
file = 'ex_im_cn_N10_50&200e-3.txt';
A = importdata(file, ',');

dt = log10(5e-3.*[1, 1./5, 1./25, 1./125, 1./625, 1./3125, 1./15625]);

figure();
hold on;
plot(dt,A(:,5),'r-')
plot(dt, (A(:,6)),'b-');
plot(dt, (A(:,7)),'c-');

ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error as a function of $\Delta t$ at $t=20\times10^{-2}$ N=10','Interpreter','LateX','Fontsize',18);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',15);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter','latex');
grid on


% file = 'ex_im_cn_an_t50e-3_dt5e-3_N10.txt';
% A = importdata(file, ',');
% figure();
% 
% x = 0:0.1:1;
% figure();
% hold on;
% plot(x, A(5,:),'r-')
% plot(x, A(6,:),'-')
% plot(x, (A(7,:)),'-');
% plot(x, (A(8,:)),'-');

%% dx = 1/100
% in x = 0.8
%% t = 50e-3
file = 'ex_im_cn_N100_50&200e-3.txt';
A = importdata(file, ',');

dt = A(:,1);
figure();
hold on;
plot(dt,A(:,2),'r-')
plot(dt, (A(:,3)),'b-');
plot(dt, (A(:,4)),'c-');

ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error as a function of $\Delta t$ at $t=5\times10^{-2}$ N=100','Interpreter','LateX','Fontsize',18);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',15);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter','latex');
grid on


% file = 'ex_im_cn_an_t50e-3_dt1e-5_N100.txt';
% A = importdata(file, ',');
% 
% x = 0:0.01:1;
% figure();
% hold on;
% plot(x, A(1,:),'r-')
% plot(x, A(2,:),'k-')
% plot(x, (A(3,:)),'b-');
% plot(x, (A(4,:)),'c-');


%% t=200e-3
hold on;
plot(dt,A(:,5),'r-')
plot(dt, (A(:,6)),'b-');
plot(dt, (A(:,7)),'c-');

ll=legend(gca,'show','explicit','implicit','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error as a function of $\Delta t$ at $t=20\times10^{-2}$ N=100','Interpreter','LateX','Fontsize',18);
xlabel('$\log_{10}(\Delta t)$','interpreter','latex','fontsize',15);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter','latex');
grid on



% file = 'ex_im_cn_an_t50e-3_dt1e-5_N100.txt';
% A = importdata(file, ',');
% 
% x = 0:0.01:1;
% figure();
% hold on;
% plot(x, A(5,:),'r-')
% plot(x, A(6,:),'k-')
% plot(x, (A(7,:)),'b-');
% plot(x, (A(8,:)),'c-');
% 











%% explicit method failure at t fixed = 0.6032
% dt = 0.501 and N=100 
file = '1204xdt0_501-N100-explicit.txt'
E = importdata(file, ',');

hold on
plot(0:0.01:1,E(2,:));
plot(0:0.01:1,E(1,:));
ylabel('u(x,t)');
xlabel('x')
ll = legend('analytical solution', 'explicit method');
set(ll);




%% explicit and implicit methods error as a function of dt
% at x = 0.9
% N = 100
% x = log10(dt)         y=log10(relative error)
file = 'error_imp_cn_t2e-3_1e-2_N100_x0_9.txt';
A = importdata(file, ',');
%% t1 = 2e-3
hold on;
plot(A(:,1),A(:,2),'r-')
plot((A(:,1)), (A(:,3)),'b-');
xlabel('log10(dt)');
ylabel('log10(relative error)')
ll = legend('implicit', 'crancknicolson');
set(ll);
%% t2 = 1e-2
hold on;
plot((A(:,1)), (A(:,4)),'r-.');
plot((A(:,1)), (A(:,5)),'b-.');
xlabel('log10(dt)');
ylabel('log10(relative error)')
ll = legend('implicit', 'crancknicolson');
set(ll);



%% explicit and implicit methods error as function of dx
% at x = 0.9
file = 'speranza_10.txt';
A = importdata(file, ',');
%% t1 = 2e-3
hold on;
plot(A(:,1),A(:,2),'r-','Linestyle','none','Marker','x')
plot((A(:,1)), (A(:,3)),'b-','Linestyle','none','Marker','x');
xlabel('log10(dx)');
ylabel('log10(relative error)')
ll = legend('implicit', 'crancknicolson');
set(ll);
%% t2  = 1e-2
figure();
hold on;
plot((A(:,1)), (A(:,4)),'r-.');
plot((A(:,1)), (A(:,5)),'b-.');
xlabel('dx');
ylabel('log(relative error)')
ll = legend('implicit', 'crancknicolson');
set(ll);



%% 2d solution at x,y = 0.9 and t = 0.1
% as a function of dt and dx
% log10 dt
ldt = -4:1:-2;
% log dx 
ldx = -3:1:-1;

file = 'asticazzi2.txt';
B = importdata(file, ',');
hold on;
plot(1e-3:1e-3:4.999, B(1,:));
plot(1e-3:1e-3:4.999, B(2,:));

ll=legend(gca,'show','N=100','N=10','crank nicolson');
set(ll,'interpreter','latex','fontsize',14,'location','northwest');
grid on;
title('Relative error as a function of time with $\Delta t=10^{-3}$ at x=y=0.5','Interpreter','LateX','Fontsize',18);
xlabel('$t$','interpreter','latex','fontsize',15);
ylabel('$\log_{10}(\epsilon)$','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter','latex');
grid on

%
%surf(B)
%colormap('parula(5)')
%%
figure();
xx = 0.01:0.01:0.99;
for ii = 0:1:max(size(B));
    
    q = 99*ii;
    qq = (ii+1)*99;
    surf(xx',xx',B((q+1):(qq),1:99));
    hold on
    %zlim([0 1])
    pause(0.0001)
    clf
end



%%
file = 'prova';
B = importdata(file, ',');
for ii = 1:1:3000000
    
    hold on
    plot(0:0.01:1, (B(ii,:)),'g')

    
    %axis([0 1 0 3e-3]);
    %plot(100:10:1000,A(1,ii));
    pause(0.1)
    clf
end









%% 2d solution difference between the two solution
% as a function of spatial coordinate
file = 'test.txt';
B = importdata(file, ',');
%%
plot(0:0.001:1, B(2,:),'rd');
hold on;
plot(0:0.001:1, B(1,:),'bx');
%%
file = 'ex.txt';
A = importdata(file, ',');


for ii = 1:10:3000000
    
    hold on
    plot(0:0.01:1, A(ii,:))
    %axis([0 1 -1 0]);
    %plot(100:10:1000,A(1,ii));
    pause(0.01)
    clf
end

%% analytical solution
file = 'stab.txt';
B = importdata(file, ',');
figure();

for ii = 1:10:max(size((B)))
    
    hold on
    plot(0:0.01:1, B(ii,:))
    %axis([0 1 0 1]);
    %plot(0:0.01:1,A(ii,:));
    pause(0.001)
    clf
end

%% 2d at x =0.5 and t = 0.01
% dt = 1e-3
% as a function of N=10->1000
file = 'errordifference2d.txt';
B = importdata(file, ',');
dt = 4:-1:1;
dx = 3:-1:1;

surf(dt,dx,B)

%%
%plot(0:0.01:1,A(:,100))
figure();
%zlim('manual',[0 1])
xx = 0:0.01:1;
for ii = 0:1:max(size(B));
    
    q = 101*ii;
    qq = (ii+1)*101;
    surf(xx',xx',B((q+1):(qq),1:101));
    hold on
    pause(0.1)
    clf
end
%%
figure();
xx = 0.01:0.01:0.99;
for ii = 1:10:max(size(B));
    
    q = 99*ii;
    qq = (ii+1)*99;
    surf(xx,xx,B((q+1):(qq),1:99));
    hold on
    %zlim([0 1])
    pause(0.0001)
    clf
end
