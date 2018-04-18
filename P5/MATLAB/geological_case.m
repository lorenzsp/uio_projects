%% Geological application of heat equation

%% Declaration of constants
from_ms_Gy=(365*24*3600*1e9);
k=2.5;
cp=1000;
rho=3.5e3;
cost_time=(rho*cp)*1e6/(k*from_ms_Gy);
time_resol=100000;
final_time=1/cost_time;
initial_time=0;
step_time=(final_time-initial_time)/(time_resol);
%% Save images for video
close all;
clear s;
f=0;
for j=0:1:size(u)
   close all
   q=121*(j);
   q1=q+20;
   q2=q1+20;
   qq=(j+1)*121;
   hold on;
plot(-y(q+1:qq),u(q+1:qq));
pause(1);
t=j*step_time*cost_time;
g=sprintf('title(''Heat diffusion in crust at $t=$%d$ Gy$'',''Interpreter'',''LateX'',''Fontsize'',18)', t);
eval(g);
xlabel('$y\ [Km]$','interpreter','latex','fontsize',15);
ylabel('$T\ [^oC]$','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter','latex');

ll=legend(gca,'show','data');
set(ll,'interpreter','latex','fontsize',14,'location','northeast');
grid on;
set(gca, 'FontName', 'Times');
set(gca,'FontSize',14);
if(f<10)
   s=sprintf('print(fig,''image2_00%d'',''-djpeg'')', f);
elseif (f<100) 
   s=sprintf('print(fig,''image2_0%d'',''-djpeg'')', f);
elseif (f<1000)
   s=sprintf('print(fig,''image2_%d'',''-djpeg'')', f);
end
   f=f+1;
   eval(s);
   pause(1);
   clf
end

%% Plots
hold on
clear s;
f=0;
fig=figure();
j1=0; %100
q=121*(j1);
qq=(j1+1)*121;
plot(-y2(q+1:qq),u2(q+1:qq));
hold on
% t=j*step_time*cost_time;
% j2=1000;
% q=121*(j2);
% qq=(j2+1)*121;
% plot(-y2(q+1:qq),u2(q+1:qq));
j4=30000;
q=121*(j4);
qq=(j4+1)*121;
plot(-y2(q+1:qq),u2(q+1:qq));
j3=100000;
q=121*(j3);
qq=(j3+1)*121;

plot(-y2(q+1:qq),u2(q+1:qq));
t=j*step_time*cost_time;

title('Heat diffusion in crust without radioactive decay','Interpreter','LateX','Fontsize',18);
xlabel('$y\ [Km]$','interpreter','latex','fontsize',15);
ylabel('$T\ [^oC]$','interpreter','latex','fontsize',15);
set(gca,'TickLabelInterpreter','latex');
t1=j1*step_time*cost_time
%t2=j2*step_time*cost_time
t4=j4*step_time*cost_time
t3=j3*step_time*cost_time %'Data $t=0.01\ Gy$' 1\times 10^{-3},
ll=legend(gca,'show','Data $t=0\ Gy$','Data $t=0.3\ Gy$','Data $t=1\ Gy$');
set(ll,'interpreter','latex','fontsize',14,'location','northeast');
grid on;
set(gca, 'FontName', 'Times');
set(gca,'FontSize',14);