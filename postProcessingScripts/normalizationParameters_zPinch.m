%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%   script to compute normalized values for MHD simulation of m=0
%%%   mode on a zPinch equilibrium
%%%
%%%   Using SI units
%%%
%%%   Specified scales: N0 [1/m^3], T0 [eV], r0 [m]
%%%
%%%   Derived scales: 
%%%   P0 = 2.0*N0*T0*qe  [J/m^3]
%%%   B0 = sqrt(P0*mu0)  [T]
%%%   J0 = B0/r0/mu0     [A/m^2]
%%%   rho0 = Mi*N0;      [kg/m^3]
%%%   U0 = sqrt(P0/rho0) [m/s]
%%%   E0 = U0*B0         [V/m]
%%%   t0 = r0/U0         [s]
%%%
%%%   I0 =
%%%   B0 = I0*mu0/r0, P0 = B0^2/mu0
%%%   
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


%%%   specify spatial scale, gas temp/pressure, atomic mass number, 
%%%   and current scale
%
aMn = 1.007;      % atomic mass number
Tg  = 300;    % ambient gas temperature [K]
Pg  = 1;      % ambient gas pressure    [Torr]
N0 = 1.0e24; %2*2.6868e25*Pg/760*273/Tg;  % total density [1/m^3]
r0  = 5.0e-3/3; % spatial scale [m]
T0  = 100.0;    % initial plasma temperature [eV]
%
aMn = 1.007;      % fuze deuterium
r0 = 0.091e-2;  % fuze a=0.09 cm,
N0 = 4.25e24;   % fuze n0 = 4e18/cc
T0 = 1270;      % fuze temperature [eV]

% aMn = 2.0;     % zap deuterium
% r0 = 1.0e-2;   % zap a=1.0 cm, r0 = 3*a
% N0 = 1.0e22;   % zap n0 = 1e16/cc
% T0 = 100;      % zap [eV]

% aMn = 2.0;     % reactor deuterium
% r0 = 5.0e-5;   % reactor 
% N0 = 3.0e27;   % reactor
% T0 = 4e4;      % reactor [eV]



%%%   fundamental constants
%
kB  = 1.3807e-23;      % Boltzmann constant [J/K]
qe  = 1.602176487e-19; % electron charge [C]
amu = 1.6605402e-27;   % atomic mass unit [kg]
me  = 9.1094e-31;      % electron mass [kg]
mu0 = 4*pi*1e-7;       % permeability of free space [H/m]
epsilon0 = 8.8542e-12; % permittivity of free space [F/m]
cvac = 2.9979e8;       % speed of light [m/s]


%%%   calculate characteristic values
%
P0 = 2.0*N0*T0*qe;            % pressure [J/m^3]
B0 = sqrt(mu0*P0);            % magnetic field [T]
%B0 = I0*mu0/r0;               % magnetic field [T]
%P0 = B0^2/mu0;                % pressure [kg/s^2/m]
J0 = B0/r0/mu0;               % current density [A/cm^2]
Mi = aMn*amu;                % particle mass [kg]
rho0 = Mi*N0;                 % density [kg/m^3]
U0 = sqrt(P0/rho0);           % velocity [m/s]
E0 = U0*B0;                   % electric field [V/m]
t0 = r0/U0;                   % time [s]


%%%   calculate fundamental values
%
wci = qe*B0/Mi; % ion cyclotron frequency [Hz]
wce = qe*B0/me; % electron cyclotron frequency [Hz]
VTe = 4.19e7*sqrt(T0)/100;
VTi = VTe*sqrt(me/Mi);
wpe = 5.64e4*sqrt(N0*1e-6);   % electron plasma freq [rad/s]
wpi = wpe*sqrt(me/Mi);        % ion plasma freq [rad/s]


% %%%   set parameters for calculation of time-dependent magnetic field
% %%%   boundary: B|x=0 = mu*I(t)/dy, I(t) = I0*t/100ns
% %
% I0  = 1e5;         % current magnitude [A]
% dy  = 2*pi*1.5e-2; % in-plane thickness [m]

% %%%   calculate scale for magnetic field bounadry condition
% %%%   Bnorm|x=0 = Bx0*tnorm/tx0;
% %
% Bx0 = mu0*I0/dy/B0;
% tx0 = 200e-9/t0;


%%%   calculate dimensionless parameters
%
delta = U0^2/cvac^2;
eta = 1.03e-4/10/T0^1.5;  % plasma resistivity [Ohm-m]
eta0  = r0^2*mu0/t0;
etanorm = eta/eta0;
rhoi = VTi/wci;                % ion gyro radius [m]


%%%   calculate collision time (Braginskii pg 215 (1965))
%
Clog = 24-log(sqrt(N0*1e-6)/T0);
Clogii = 23-log(sqrt(2*N0*1e-6)/T0^1.5);
tau_e = 3.5e4*T0^1.5/(Clog/10)/(N0*1e-6); % electron collison time [s]
tau_i = 3.0e6*T0^1.5/(Clogii/10)/(N0*1e-6)*sqrt(aMn/2); % ion collison time [s]

nuT = me/Mi/tau_e; % thermalization rate [Hz]

%%%   calculate ion and electron inertial length scales (skin depth)
%
Le = cvac/wpe;                % electron inertial scale [m]
Li = cvac/wpi;                % ion inertial scale [m]


%%%   calculate gyro-Bohm radius
%
Cs = sqrt(P0/rho0);
rhos = Cs/wci; % [m]

display(delta);
display(etanorm);


