%%%%%%%%%%%%%%%%%%%%
%%%
%%%    d(P+B^2/2)/dr + B^2/r  = 0
%%%    or 
%%%    d(rP)/r - P/r + d(rB)/dr/r*B = 0;
%%%
%%%    d(r*(P+B^2/2))/dr - (P-B^2/2) = 0
%%%
%%%    calculate stability region from Kadomtsev 1960 Eqs 29 and 31
%%%    for Bennett equilibrium profiles of zpinch
%%%
%%%    P = a^4/(a^2+r^2)^2
%%%    B = sqrt(2)*a*r/(a^2+r^2)
%%%
%%%    gamma = 5/3, beta = 2*P/B^2
%%%
%%%    ideal sausage: w^2 ~ dlnP/dlnr + 4*gamma/(2+gamma*beta)
%%%
%%%    drift-ideal entropy (if sausage satisifed)
%%%    dlnT/dlnr + 1 + (7/10 + beta/4)*dlnP/dlnr > 0:
%%%    Note that 7/10 = (2*gamma-1)/(2*gamma), gamma = 5/3
%%%
%%%
%%%    r is normalized by some r0, 
%%%    U is normalized by sqrt(P0/(MiN0)) = sqrt(2*T0/Mi) = Cs0
%%%    omega is normalized by omega0=Cs0/r0
%%%
%%%%%%%%%%%%%%%%%%%%
clear all;


%%%   specify input parameters
%
gamma = 5/3; %2.1; %5/3; %2.1;  % adiabatic coefficient
ka = 30;      % normalized wavenumber
Li = 1.5e-6;  % omega0/Omegai0 = ion inertial length / r0
Li = 0.012;  % omega0/Omegai0 = ion inertial length / r0

nuT = 0; %200;    % normalized thermalization rate


%%%   set radial grid
%
nr = 400;
R  = 3;
a  = 1;
dr = R/nr;


%%%   calculate Bennett profiles, qP, and plasma beta
%
rcc = dr/2:dr:R-dr/2;
x = rcc/a;
P = 1./(1+x.^2).^2;
B = sqrt(2)*x./(1+x.^2);
P(end) = P(end-1);
B(end) = B(end-1)*rcc(end-1)/rcc(end);
qP = -4*x.^2./(1+x.^2); % dlnP/dlnr
qT = 0*x;               % dlnT/dlnr
beta = 1./x.^2;         % plasma beta, 2*P/B^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%          compute modes in limit where k*rhos = kR*Li << 1
%%%          Kadomtsev 1960 Eqs 29 and 30
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   compute growth rate from ideal MHD (neglect magnetosonic waves)
%
omega2 = 2.0./rcc.^2.*(qP + 4*gamma./(2 + gamma*beta));
omega = sqrt(omega2);


%%%   compute growth rate from entropy mode (Kadomtsev)
%
nuTstar = 2*nuT*rcc.*B/(ka*Li);
a0 = (2*rcc.*B/(ka*Li)).^2;
b0 = -2*rcc.*B/(ka*Li).*(2*qP - 1i*nuTstar);
RHS = gamma^2*(4 + beta.*qP).^2./(4*gamma + qP.*(2 + gamma*beta)).* ...
      (1 + ((2*gamma-1)./(2*gamma)+beta/4).*qP - qT);
c0 = qP.*(qP -1i*nuTstar) - RHS;

omega_entropy = zeros(length(rcc),2);
for i=1:length(rcc)
    thisp = [a0(i) b0(i) c0(i)];
    omega_entropy(i,:) = roots(thisp);
end
real_entropy = qP./rcc./B/2*ka*Li;
gamma_entropy = sqrt(-RHS)./rcc./B/2*ka*Li;
%omega2_entropy = 1+((2*gamma-1)/(2*gamma)+beta/4).*qP;
%figure(4); hold on; plot(x,omega2_entropy); box on;


%%%   plot equilibrium profiles
%
close(figure(11));
f11=figure(11); set(f11,'position',[1400 900 560 420]);
plot(rcc,P,rcc,B,rcc,beta); axis([0 3 0 1]);
xlabel('r/r_0'); title('normalized equilibrium profiles');
legend('pressure','magnetic field','plasma \beta=2P/B^2');
axis('square'); grid on;


%%%   plot growth rate for ideal modes
%
f2=figure(2); set(f2,'position',[1400 400 1100 420]);
subplot(1,2,1);
hold on; plot(rcc,real(omega),'displayName','ideal'); box on; 
hold on; plot(rcc,real(-omega),'displayName','ideal'); box on; 
title('real roots'); xlabel('r/a'); ylabel('\omega/\omega_0');
%lg21 = legend('show'); set(lg21,'location','best');
%
subplot(1,2,2);
hold on; plot(rcc,imag(omega),'displayName','ideal'); box on; 
hold on; plot(rcc,imag(-omega),'displayName','ideal'); box on; 
title('imaginary roots'); xlabel('r/a'); ylabel('\gamma/\omega_0');
%lg22 = legend('show'); set(lg22,'location','best');


%%%   plot growth rate for entropy modes
%
%f3=figure(3); set(f3,'position',[1400 100 1100 420]);
%
subplot(1,2,1);
%hold on; plot(rcc,real_entropy,'displayName','entropy'); box on; 
hold on; plot(rcc,real(omega_entropy(:,1)),'displayName','entropy'); box on; 
hold on; plot(rcc,real(omega_entropy(:,2)),'displayName','entropy'); box on; 
title('real roots'); xlabel('r/a'); ylabel('\omega/\omega_0');
lg21 = legend('show'); set(lg21,'location','best');
%
subplot(1,2,2);
%hold on; plot(rcc,gamma_entropy,'displayName','entropy'); box on; 
hold on; plot(rcc,imag(omega_entropy(:,1)),'displayName','entropy'); box on; 
hold on; plot(rcc,imag(omega_entropy(:,2)),'displayName','entropy'); box on; 
title('imaginary roots'); xlabel('r/a'); ylabel('\gamma/\omega_0');
lg22 = legend('show'); set(lg22,'location','best');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%          compute modes from full dispersion (sans thermalization)
%%%          Kadomtsev 1960 Eqs 26-28
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

krhoi = ka*Li./B/sqrt(2);
krhos = krhoi*sqrt(2);

c4 = (gamma-1)/gamma-(1+beta/2);
c3 = 2*qP.*(1+beta/2-(gamma-1)/gamma) - 4;
c2 = (gamma-1)/gamma*qP.*(qP-gamma*(2+beta/2.*qP)-4./krhoi.^2) ... 
   + 4./krhoi.^2.*(2+(1+beta/2).*qP) + 8*qP ...
   - (1+beta/2).*(qP.^2-gamma*(2+beta/2.*qP).*(2+((gamma-1)/gamma+beta/2).*qP));
c1 = (gamma-1)/gamma*8*qP.^2./(krhoi).^2 - 8*qP.*(2+(1+beta/2).*qP)./krhoi.^2 ...
   - 4*qP.^2 + 4*gamma*(2+beta/2.*qP).*(2+((gamma-1)/gamma+beta/2).*qP);
c0 = -4*qP.^2./krhoi.^2.*(gamma-1)/gamma.*(qP-gamma*(2+beta/2.*qP)) ...
   + 4./krhoi.^2.*(2+(1+beta/2).*qP).*(qP.^2-gamma*(2+beta/2.*qP).*(2+((gamma-1)/gamma+beta/2).*qP));

% c4 = c4.*(2*rcc.*B/(kR*Li)).^4;
% c3 = c3.*(2*rcc.*B/(kR*Li)).^3;
% c2 = c2.*(2*rcc.*B/(kR*Li)).^2;
% c1 = c1.*(2*rcc.*B/(kR*Li)).^1;


omega_total = zeros(length(rcc),4);
for i=1:length(rcc)
    thisp = [c4(i) c3(i) c2(i) c1(i) c0(i)];
    omega_total(i,:) = roots(thisp);
    omega_total(i,:) = omega_total(i,:)*ka*Li/(2*rcc(i)*B(i));
end



%%%   plot growth from general 4th order solution
%
f3=figure(3); set(f3,'position',[1400 400 1100 420]);
%
subplot(1,2,1);
hold on; plot(rcc,real(omega_total),'displayName','general'); box on; 
title('real roots'); xlabel('r/a'); ylabel('\omega/\omega_0');
lg31 = legend('show'); set(lg21,'location','best');
%
subplot(1,2,2);
hold on; plot(rcc,imag(omega_total),'displayName','general'); box on; 
title('imaginary roots'); xlabel('r/a'); ylabel('\gamma/\omega_0');
lg32 = legend('show'); set(lg22,'location','best');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%          compute modes from MHD with full Ohms law, but
%%%          only standard entropy density equation
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c3 = (gamma-1)/gamma - (1+beta/2);
c2 = -(gamma-1)/gamma*qP - 4 + qT;
c1 = 4./krhoi.^2.*(2 + (1+beta/2).*qP - (gamma-1)/gamma*qP);
c0 = 4*qP./krhoi.^2.*((gamma-1)/gamma*qP - qT);



omega_FullOhms = zeros(length(rcc),3);
for i=1:length(rcc)
    thisp = [c3(i) c2(i) c1(i) c0(i)];
    omega_FullOhms(i,:) = roots(thisp);
    omega_FullOhms(i,:) = omega_FullOhms(i,:)*ka*Li/(2*rcc(i)*B(i));
end


%%%   plot growth from Full Ohms law
%
f4=figure(4); set(f4,'position',[1400 400 1100 420]);
%
subplot(1,2,1);
hold on; plot(rcc,real(omega_FullOhms),'displayName','Hall MHD'); box on; 
title('real roots'); xlabel('r/a'); ylabel('\omega/\omega_0');
lg31 = legend('show'); set(lg21,'location','best');
%
subplot(1,2,2);
hold on; plot(rcc,imag(omega_FullOhms),'displayName','Hall MHD'); box on; 
title('imaginary roots'); xlabel('r/a'); ylabel('\gamma/\omega_0');
lg32 = legend('show'); set(lg22,'location','best');
 


%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   plot Spitzer factors for diamagnetic heat flux
%%%
%%%



%%%   characteristic values
%
aMn = 2.0;      % fuze deuterium
r0 = 0.091e-2;   % fuze a=0.09 cm,
N0 = 4.25e24;    % fuze n0 = 4e18/cc
T0 = 1270;      % fuze [eV]


%%%   physical constants and mass calculation
%
qe  = 1.602176487e-19; % electron charge [C]
amu = 1.6605402e-27;   % atomic mass unit [kg]
me  = 9.1094e-31;      % electron mass [kg]
Mi  = aMn*amu;         % ion mass [kg]
mu0 = 4*pi*1e-7;       % permeability of free space [H/m]
%
P0 = 2.0*N0*T0*qe;            % pressure [J/m^3]
B0 = sqrt(mu0*P0);            % magnetic field [T]


%%%   calculate collision time (Braginskii pg 215 (1965))
%
Clog = 24-log(sqrt(N0*1e-6)/T0);
Clogii = 23-log(sqrt(2*N0*1e-6)/T0^1.5);
tau_e = 3.5e4*T0^1.5/(Clog/10)/(N0*1e-6); % electron collison time [s]
tau_i = 3.0e6*T0^1.5/(Clogii/10)/(N0*1e-6)*sqrt(aMn/2); % ion collison time [s]


wci = qe*B0/Mi; % ion cyclotron frequency [Hz]
wce = qe*B0/me; % electron cyclotron frequency [Hz]


%%%   calculate wctau factors for co-perp heat fluxes
%
Xele = wce*tau_e*B;
Xion = wci*tau_i*B;

Dele = Xele.^4 + 14.79*Xele.^2 + 3.7703;
Dion = Xion.^4 +  2.70*Xion.^2 + 0.6770;

coPerpFactor_ele = Xele.^2.*(5/2*Xele.^2 + 21.67)./Dele;
coPerpFactor_ion = Xion.^2.*(5/2*Xion.^2 +  4.65)./Dion;
coPerpFactor_fudge = (100*B).^2./((100*B).^2+3.7703); % simulation fudge

f11=figure(11); 
plot(rcc,coPerpFactor_ele/(5/2),'displayName','ele');
hold on; plot(rcc,coPerpFactor_ion/(5/2),'displayName','ion');
hold on; plot(rcc,coPerpFactor_fudge,'displayName','fudge');
xlabel('r/r_0'); title('co-perp factors');
lg11=legend('show');
% 






