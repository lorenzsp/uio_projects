%% kronig-penney model
% fixed constant of parametrer of the lattice a=1
% energy
E =0:0.01:20;
% hbar =1
m =1;
%strength of the potential
b = 1;
k = sqrt(2 * m *E);
f = cos(k) + b.*sin(k)./k;

figure();
plot(k,f,'.');
grid on;

count=1;
for t=1:length(k)
    if (f(t)<1) && (f(t)>-1)
        w(t)=acos(f(t));
        ew(t) = E(t);
    end
    
end
% composite function
for t=1:length(w)
    if t>495
        w(t)=6.28-w(t);
    end

    if w(t)>6.26
        w(t)=2*6.28-w(t);
    end
end

figure();
hold on;
plot(w,ew,'.');
%set(gca,'xticklabel',[])
%set(gca,'yticklabel',[])
grid on
plot_f('Energy as a function of k in the 1st BZ','k','E',14)
%xlim([0 3.15]);
%% 1st derivative
for t=1:length(ew)-1
   d_ew(t) = (ew(t+1)-ew(t)) ./(w(t+1) - w(t)); 
end
figure();
plot(w(1:end-1),d_ew,'.')
grid on
plot_f('Derivative of Energy as a function of k in the 1st BZ','k','dE/dk',14)
xlim([0 3.1]);
% 2nd derivative
for t=1:length(d_ew)-1
   dd_ew(t) = (d_ew(t+1)-d_ew(t)) ./(w(t+1) - w(t)); 
end
figure();
plot(w(1:end-2),dd_ew,'.')
grid on
plot_f('Second derivative of Energy as a function of k in the 1st BZ','k','d^2E/dk^2',14)
xlim([0 3.1]);

%% module 4 ex 3
E_0 = 1.42; % eV
m1 = 0.067; % m_e mass of the electron 0.511 MeV 
m2 = 0.082;
m3 = 0.45;

k=-0.4:0.001:0.4;
% energy = hbar^2 c^2 /(2 m c^2) k^2 = 197^2 MeV F^2 
% = 7.5947e-20 eV m = 38.809 10^-1 Angstrom^2 k^2 =3.8809 A^2 k^2
E1 = E_0 + k.^2 * 7.5947/(2*m1);
E2 = - k.^2 * 7.5947/(2*m2);
E3 = - k.^2 * 7.5947/(2*m3);

figure();
hold on;
plot(k,E1);
plot(k,E2);
plot(k,E3);
plot_f('Energy dispersion relations of GaAs', ['k ' '[' char(197) '^{-1}]'] , 'E [eV]',14)
s= {['m_1 = ',num2str(m1),' m_e'],...
    ['m_2 = ',num2str(m2),' m_e'],...
    ['m_3 = ',num2str(m3),' m_e']};
legend_f(s);


%% ex 5
%
kb = 8.617e-5; %eV/K 
hbar = 6.582e-16; %eV s
E_g = 1.12; %eV
E_d = 0.045;
m_e = 1.08 *0.511; %eV
m_h = 0.811*0.511;

T=5:1:1.1455e+03;
% cm^-3
N_d = 1e17;
N_c = 2.9e19; 
n_i = 2*N_d * ...
    (1 + ...
    sqrt( 1 + 4* (N_d / N_c) .* exp(E_d./(kb *T) )  ) ...
    ).^(-1);


figure();
plot(100./T, log10(n_i),'b');
hold on;

T1=1.1455e+03:1:0.5e4;
n_h = N_c*exp(-E_g./(2*kb*T1));
plot(100./T1, log10(n_h),'b')

plot_f('Electron concentration n as a function of temperature T','100/T [K^{-1}]','log_{10}(n)',14)
%xlim([0 20e-3]);
ylim([16 18]);


%% energy plot
T=0.0001:0.1:500;
n_i = 2*N_d * ...
    (1 + ...
    sqrt( 1 + 4* (N_d / N_c) .* exp(E_d./(kb *T) )  ) ...
    ).^(-1);
E_f = kb * T .*log(N_d./n_i);

plot(100./T, E_f,'r')
hold on;
plot(100.*[0 2e-3],[4.1615e-04 4.1615e-04],'r')
plot_f('Fermi energy E_f as a function of temperature T','100/T [K^{-1}]','E_f - E_i [eV]',14)


