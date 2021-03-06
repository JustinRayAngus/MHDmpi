%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%   script to compute normalized values for MHD simulation
%%%
%%%   Using SI units
%%%
%%%   Specified scales: N0 [1/m^3], T0 [eV], r0 [m]
%%%
%%%   Derived scales: 
%%%   P0 = N0*T0*qe      [J/m^3]
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
aMn = 2;      % atomic mass number
Tg  = 300;    % ambient gas temperature [K]
Pg  = 3;      % ambient gas pressure    [Torr]
N0 =  2*2.6868e25*Pg/760*273/Tg;  % total density [1/m^3]
T0  = 1.0;    % initial plasma temperature [eV]


%%%   set parameters for calculation of time-dependent magnetic field
%%%   boundary: B|x=0 = mu*I(t)/dy, I(t) = I0*t/100ns
%
I0  = 1e5;         % current magnitude [A]
dy  = 2*pi*1.5e-2; % in-plane thickness [m]
r0  = 1.0e-2; % spatial scale [m]

I0  = 1e6;       % current magnitude [A]
r0  = 3.0e-2;    % spatial scale [m]
dy  = 2*pi*r0;   % in-plane thickness [m]


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
J0 = B0/r0/mu0;               % current density [A/m^2]
Mi = aMn*amu;                % particle mass [kg]
rho0 = Mi*N0;                 % density [kg/m^3]
U0 = sqrt(P0/rho0);           % velocity [m/s]
E0 = U0*B0;                   % electric field [V/m]
t0 = r0/U0;                   % time [s]


%%%   calculate scale for magnetic field bounadry condition
%%%   Bnorm|x=0 = Bx0*tnorm/tx0;
%
Bx0 = mu0*I0/dy/B0;
tx0 = 200e-9/t0;


%%%   calculate ion and electron inertial length scales (skin depth)
%
wpe = 5.64e4*sqrt(N0*1e-6);   % electron plasma freq [rad/s]
wpi = wpe*sqrt(me/Mi);        % ion plasma freq [rad/s]
Le = cvac/wpe;                % electron inertial scale [m]
Li = cvac/wpi;                % ion inertial scale [m]
eta = 1.03e-4/10/T0^1.5;  % plasma resistivity [Ohm-m]


%%%   calculate dimensionless parameters
%
lambdae0sq = (Le/r0)^2;
delta0 = U0^2/cvac^2;
eta0 = eta/(r0^2*mu0/t0);



