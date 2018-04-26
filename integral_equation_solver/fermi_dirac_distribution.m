%% plot fermi dirca distribution
% mu = eps_f chemical pot equal fermi energy
r_T = [0.01, 0.1, 0.5, 1, 1.2]; % ratio between T and T_f
% energy interval E = eps/eps_f
E = 0:0.001:5;
% Fermi dirac distr
%f = 1./(1 + exp(r_T (E-1))) 
figure();
for i=1:length(r_T)
    hold on;
    f = 1./(1 + exp( (E-1) ./ r_T(1,i)));
    pp = plot(E, f,'LineWidth',2);
    
end
plot_f('\fontsize{16}Fermi-Dirac Distribution with \mu=\epsilon_F', ...
    '\fontsize{16}\epsilon/\epsilon_F', ...
    '\fontsize{16}f', 15);

s= {['T_1 = ',num2str(r_T(1,1)),' T_F'],...
    ['T_2 = ',num2str(r_T(1,2)),' T_F'],...
    ['T_3 = ',num2str(r_T(1,3)),' T_F'],...
    ['T_4 = ',num2str(r_T(1,4)),' T_F'],...
    ['T_5 = ',num2str(r_T(1,5)),' T_F']};
legend_f(s)
%% plot FDD
file = 'temperature_mu';
m = importdata(file,',');
plot(m(:,1),m(:,2),'-');
plot_f('\fontsize{16}Chemical potential', ...
    '\fontsize{16}T/T_F', ...
    '\fontsize{16}\mu/\epsilon_F', 15)
%% 
r_T = [0.01, 0.1, 0.5, 1, 1.2];% ratio between T and T_f
mu = [0.99812, 0.987, 0.743, -0.022, -0.431];
% energy interval E = eps/eps_f
E = 0:0.001:5;
% Fermi dirac distr
%f = 1./(1 + exp(r_T (E-1))) 
figure();
for i=1:length(r_T)
    hold on;
    f = 1./(1 + exp( (E-mu(1,i)) ./ r_T(1,i)));
    pp = plot(E, f,'LineWidth',2);
    
end
plot_f('\fontsize{16}Fermi-Dirac Distribution', ...
    '\fontsize{16}\epsilon/\epsilon_F', ...
    '\fontsize{16}f', 15);

s= {['T_1 = ',num2str(r_T(1,1)),' T_F'],...
    ['T_2 = ',num2str(r_T(1,2)),' T_F'],...
    ['T_3 = ',num2str(r_T(1,3)),' T_F'],...
    ['T_4 = ',num2str(r_T(1,4)),' T_F'],...
    ['T_5 = ',num2str(r_T(1,5)),' T_F']};
legend_f(s)




