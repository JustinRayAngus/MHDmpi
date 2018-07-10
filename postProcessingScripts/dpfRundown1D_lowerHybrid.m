%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   calculate lower hybrid resistivity from 1D dpf rundown sims
%%%
%%%   resistivity: eta = 4pi*nu/wpe^2 [Hz]
%%%
%%%   classical: nu = 1/taue, taue = 3.44e5*Te^1.5/N
%%%   lower-hybrid: nu = (Vd/VTi)^2*wlh, wlh = sqrt(wci*wce)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%%   set default font and lines
%
set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'defaultaxeslinewidth',1.5);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultaxesfontweight','bold');


numProcs = 4;
filePath = '../physicsMods/dpfRundown1D/'; TwoTempVersion=0;
filePath = '../physicsMods/dpfRundown1D/2TempVersion/'; TwoTempVersion=1;
filePath = '../physicsMods/dpfRundown1D/2TempVersion/dataSave_1MAcyl/'; TwoTempVersion=1;

xp1 = 1;

Xcc = loadData1DVec(filePath,numProcs,'Xcc');
Xce = loadData1DVec(filePath,numProcs,'Xce');
tout = loadData1DVec(filePath,numProcs,'tout');


N  = loadData1DVec(filePath,numProcs,'N');
M  = loadData1DVec(filePath,numProcs,'M');
B  = loadData1DVec(filePath,numProcs,'B');
P  = loadData1DVec(filePath,numProcs,'P');
if(TwoTempVersion)
    Pe  = loadData1DVec(filePath,numProcs,'Pe');
    Pi  = loadData1DVec(filePath,numProcs,'Pi');
    Te  = loadData1DVec(filePath,numProcs,'Te');
    Ti  = loadData1DVec(filePath,numProcs,'Ti');
else
    T  = loadData1DVec(filePath,numProcs,'T');
    Te = T/2;
    Ti = T/2;
    Pe = P/2;
    Pi = P/2;
end
V  = loadData1DVec(filePath,numProcs,'V');
J  = loadData1DVec(filePath,numProcs,'J');
Jcc = loadData1DVec(filePath,numProcs,'Jcc');
J0  = loadData1DVec(filePath,numProcs,'J0');
Ez  = loadData1DVec(filePath,numProcs,'Ez');
eta = loadData1DVec(filePath,numProcs,'eta');
gamma0 = loadData1DVec(filePath,numProcs,'gamma0');
delta0 = loadData1DVec(filePath,numProcs,'delta0');
FluxM = loadData1DVec(filePath,numProcs,'FluxM');


%%%   get scales from output files
%
tscale = loadData1DVec(filePath,numProcs,'tscale'); % [s]
Xscale = loadData1DVec(filePath,numProcs,'Xscale'); % [m]
Vscale = loadData1DVec(filePath,numProcs,'Vscale'); % [m/s]
Nscale = loadData1DVec(filePath,numProcs,'Nscale'); % [1/m^3]
Tscale = loadData1DVec(filePath,numProcs,'Tscale'); % [eV]
Bscale = loadData1DVec(filePath,numProcs,'Bscale'); % [T]
Ezscale = loadData1DVec(filePath,numProcs,'Ezscale'); % [V/m]
Jscale = loadData1DVec(filePath,numProcs,'Jscale'); % [A/m^2]
Mi = loadData1DVec(filePath,numProcs,'Mi'); % [kg]
me = 9.1094e-31; % [kg]
qe = 1.6022e-19; % [C]

%%%   calculate characteristic frequency scales
%
wpe = 5.64e4*sqrt(N*Nscale/1e6); % electron plasma frequency [rad/s]
wpi = wpe*sqrt(me/Mi);           % ion plasma frequency [rad/s]
wce = qe*abs(B)*Bscale/me;            % electron cyclotron frequency [rad/s]
wci = qe*abs(B)*Bscale/Mi;            % ion cyclotron frequency [rad/s]
wlh0 = sqrt(wce.*wci);             % lowest order lower-hybrid frequency [rad/s]
wlh  = wlh0.*wpi./sqrt(wpi.^2 + wlh0.^2);             % lower-hybrid frequency [rad/s]
%
taue = 3.44e5/10*(Te*Tscale).^1.5./(N*Nscale/1e6); % electron collision time [s]
Vez  = Jcc*Jscale/qe./(N*Nscale)*100; % electron drift speed [cm/s]
VTe  = 4.2e7*sqrt(Te*Tscale);       % electron thermal speed [cm/s]
VTi  = VTe*sqrt(me/Mi);             % ion thermal speed [cm/s]


nuei = 1./taue;              % classical collision frequency [Hz]
nulh0 = (Vez./VTi).^2.*wlh0;    % lowest order lower-hybrid collision frequency [Hz]
nulh = (Vez./VTi).^2.*wlh;      % lower-hybrid collision frequency [Hz]


f1 = figure(4); set(f1,'position',[1300 200 1100 900]);
it = 100;
%
subplot(2,2,1); 
hold on; plot(Xcc*Xscale*100,N(:,it)*Nscale/1e6); box on;
xlabel('x [cm]'); ylabel('density [1/cm^3]'); grid on;
title('plasma density'); xlim([0 1*Xscale*100]);
%
subplot(2,2,2); 
hold on; plot(Xce*Xscale*100,J(:,it)*Jscale*1e-7); box on;
xlabel('x [cm]'); ylabel('current density [kA/cm^2]'); grid on;
title('current density'); xlim([0 1*Xscale*100]);
%
subplot(2,2,3); 
hold on; plot(Xcc*Xscale*100,Te(:,it)*Tscale,'displayName','ele'); box on;
hold on; plot(Xcc*Xscale*100,Ti(:,it)*Tscale,'displayName','ion'); 
xlabel('x [cm]'); ylabel('temperature [eV]'); grid on;
title('plasma temperature'); xlim([0 1*Xscale*100]);
lg4 = legend('show'); set(lg4,'location','best');
%
subplot(2,2,4); 
hold on; plot(Xcc*Xscale*100,nuei(:,it),'displayName','spitzer'); box on;
hold on; plot(Xcc*Xscale*100,nulh(:,it),'displayName','lower-hybrid');
hold on; plot(Xcc*Xscale*100,nulh0(:,it),'linestyle','--','displayName','lower-hybrid_0');
set(gca,'yscale','log');
xlabel('x [cm]'); ylabel('frequency [Hz]'); grid on; grid off; grid on;
title('collison frequencies'); 
xlim([0 1*Xscale*100]); ylim([1e6 1e12]);
lg4 = legend('show'); set(lg4,'location','northeast');



