%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 1D Sod shock simulations
%%%
%%%   see G.A. Sod,. J. Comp. Phys. 27, 1-31 (1978)
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
filePath = '../physicsMods/dpfRundown1D_2Temp/dataSave_1MAcyl/';
filePath = '../physicsMods/dpfRundown1D_2Temp/dataSave_1MAcar/';
filePath = '../physicsMods/dpfRundown1D_2Temp/dataSave_200kAcyl/';
filePath = '../physicsMods/dpfRundown1D_2Temp/';

plotBackIndex = 0*35; % plot time will be end-plotBackIndex
xp1 = 1;

Xcc = loadData1DVec(filePath,numProcs,'Xcc');
Xce = loadData1DVec(filePath,numProcs,'Xce');
tout = loadData1DVec(filePath,numProcs,'tout');

% fileName = ['output',num2str(i-1),'.h5'];
% thisFile = [filePath,fileName];
% procID  = hdf5read(thisFile,'procID');
% fileinfo = hdf5info(thisFile);

Btot = 0*tout;
Poynt = 0*tout;
Bflux = 0*tout;
Mflux = 0*tout;
Msource = 0*tout;
Ssource = 0*tout;
Psource = 0*tout;
Poynt_exit = 0*tout;
EEztot = 0*tout;

N  = loadData1DVec(filePath,numProcs,'N');
M  = loadData1DVec(filePath,numProcs,'M');
B  = loadData1DVec(filePath,numProcs,'B');
P  = loadData1DVec(filePath,numProcs,'P');
if(TwoTempVersion)
    Pe  = loadData1DVec(filePath,numProcs,'Pe');
    Pi  = loadData1DVec(filePath,numProcs,'Pi');
    Te  = loadData1DVec(filePath,numProcs,'Te');
    Ti  = loadData1DVec(filePath,numProcs,'Ti');
    %
    hy_cc = loadData1DVec(filePath,numProcs,'hy_cc');
    hy_ce = loadData1DVec(filePath,numProcs,'hy_ce');
else
    T  = loadData1DVec(filePath,numProcs,'T');
    Te = T/2;
    Ti = T/2;
    Pe = P/2;
    Pi = P/2;
    %
    hy_cc = 1.0 + 0.0*Xcc;
    hy_ce = 1.0 + 0.0*Xce;
end
Te = 2.0*Te; Ti = 2.0*Ti; % normalization is 2*T0
T = (Te+Ti)/2.0;
V  = loadData1DVec(filePath,numProcs,'V');
J  = loadData1DVec(filePath,numProcs,'J');
Jcc = loadData1DVec(filePath,numProcs,'Jcc');
J0  = loadData1DVec(filePath,numProcs,'J0');
Ez  = loadData1DVec(filePath,numProcs,'Ez');
eta = loadData1DVec(filePath,numProcs,'eta');
gamma0 = loadData1DVec(filePath,numProcs,'gamma0');
delta0 = loadData1DVec(filePath,numProcs,'delta0');
FluxM = loadData1DVec(filePath,numProcs,'FluxM');
FluxN = loadData1DVec(filePath,numProcs,'FluxN');

%T = P/2.0./N; % average temperature

%%%   calculate divU
%
dX = Xcc(2)-Xcc(1);
dV = dX*hy_cc;
divV = zeros(size(V));
dPdx = zeros(size(V));
dhydx = zeros(size(Xcc));
for j=2:length(Xcc)-1
    divV(j,:) = (V(j+1,:)-V(j-1,:))/2/dX;
    dPdx(j,:) = (P(j+1,:)-P(j-1,:))/2/dX;
    dhydx(j) = (hy_ce(j)-hy_ce(j-1))/dX;
end
dhydx(1) = dhydx(2);
dhydx(end) = dhydx(end-1);


%%%   calcualte Ez at cell-center
%
Ezcc = zeros(size(N));
for j=3:length(Xcc)-2
    Ezcc(j,:) = (Ez(j-1,:) + Ez(j,:))/2;
end


%figure(77); hold on; plot(Xcc,divV(:,100)); box on;

Pmag = B.^2/2;
Pmageff = Pmag + cumtrapz(Xcc,B.^2./Xcc); % for cyl only
Pram = M.*V;
Ptot = P + B.^2/2 + M.*V/2.0;
PdV = P.*divV;
VdP = V.*dPdx;
etaJ2 = eta.*Jcc.^2;
etaJ2mod = etaJ2./Te;
Edenstot = 3/2*P + N.*V.^2/2 + B.^2/2;


%%%   calculate <Te+Ti> in compressed region at stagnation
%
mu0 = 4*pi*1e-7;
qe = 1.6022e-19;
I0 = loadData1DVec(filePath,numProcs,'Iscale'); % current scale [A]
N0 = loadData1DVec(filePath,numProcs,'Nscale'); % dens scale [1/m^3]
R0 = loadData1DVec(filePath,numProcs,'Xscale'); % length scale [m]
rs = 0.13*R0; % compression radius from 1D shock
TeTi0 = (gamma0-1)/(pi*R0^2*N0)*mu0*I0^2/(4*pi)*log(R0/rs)/qe; % <Ti+Te> [eV]


%%%   calculate conservation stuff
%
JdotE  = sum(Jcc(3:end-2,:).*Ezcc(3:end-2,:).*dV(3:end-2));
Etot   = sum(Edenstot(3:end-2,:).*dV(3:end-2));
Etherm = sum(1.5*P(3:end-2,:).*dV(3:end-2));
Emean  = sum(0.5*M(3:end-2,:).*V(3:end-2,:).*dV(3:end-2));
Efield = sum(B(3:end-2,:).^2/2.*dV(3:end-2));
EEztot = delta0*sum(Ez(2:end-2,:).^2/2.0.*dV(3:end-2));
%

EntropyDensity = N.*log(Te.^1.5./N) + N.*log(Ti.^1.5./N);
Entropy  = sum(EntropyDensity(3:end-2,:).*dV(3:end-2));
Momentum = sum(M(3:end-2,:).*dV(3:end-2));
Mass     = sum(N(3:end-2,:).*dV(3:end-2));

for n=1:length(tout)

    Btot(n)   = sum(B(3:end-2,n))*dX;
    Ssource(n) = sum(etaJ2mod(3:end-2,n).*dV(3:end-2));
    Msource(n) = sum(-Jcc(3:end-2,n).*B(3:end-2,n).*dV(3:end-2) ...
                     +P(3:end-2,n).*dhydx(3:end-2)*dX);
    Psource(n) = sum(etaJ2(3:end-2,n)-0*PdV(3:end-2,n).*dV(3:end-2));
   % Psource(n) = sum(etaJ2(3:end-2,n)+VdP(3:end-2,n))*dX;
    Etot(n)    = Etot(n) + EEztot(n);
    
    Poynt(n) = Ez(end-1,n)*(B(end-2,n)+B(end-1,n))/2.0 ...
             - Ez(2,n)*(B(3,n)+B(2,n))/2.0 ...
             - 0*5/2*(P(2,n)+P(3,n))/2*(V(2,n)+V(3,n))/2 ...
             - 0*0.5*N(3,n)*V(3,n)^3;
    Bflux(n) = Ez(end-1,n) - Ez(2,n);
    Mflux(n) =  (P(2,n)+P(3,n))/2 + B(2,n)^2/2;
    Msource(n) = Msource(n) + FluxM(2,n) - FluxM(end-1,n);
  %  Psource(n) = Psource(n) - 3/2*P(2,n)*(V(2,n)+V(3,n))/2;

    Msource(n) = Msource(n) - (P(end-2,n)+P(end-1,n))/2;
    Mflux(n) = Mflux(n) - (P(end-2,n)+P(end-1,n))/2;
   
end
Cs = sqrt(gamma0*P./N);
Mach = abs(V)./Cs;


% rho2/rho1 = v1/v2 = (gamma0+1)/(gamma0-1+2/M1^2)
M1 = 1.1;
deltaN = (gamma0+1)/(gamma0-1+2/M1^2);
deltaP = (2*gamma0*M1.^2-(gamma0-1))/(gamma0+1);

Ez0 = eta.*Jcc-V.*B;


% %%%
% figure(10); hold on; plot(Xcc,Ptot(:,end),'black','displayName','total');
% hold on; plot(Xcc,B(:,end).^2/2.0,'b','displayName','magnetic');
% hold on; plot(Xcc,1.5*P(:,end),'r','displayName','thermal');
% hold on; plot(Xcc,M(:,end).*V(:,end)/2.0,'color',[.47 .67 .19],'displayName','mean');
% title('total pressure');
% box on;
% if(i==1) 
%     lg10=legend('show'); set(lg10,'location','best');
% end


plots=1;
if(plots)
f1=figure(11); 
set(f1,'position',[1030 425 1300 840]);
%set(f1,'position',[341 436 900 840]);

subplot(2,3,1);
hold on; plot(Xcc,N(:,1),'black'); box on;
hold on; plot(Xcc,N(:,round((end-plotBackIndex)/2)),'b');
%hold on; plot(Xcc,N(:,round(120)),'g');
%hold on; plot(Xcc,N(:,round(end/2)),'b');
hold on; plot(Xcc,N(:,end-plotBackIndex),'r'); grid on;
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('N');
title('mass density'); axis('square');
xlim([0 xp1]);
%
subplot(2,3,2);
hold on; plot(Xcc,V(:,1),'black'); box on;
hold on; plot(Xcc,V(:,round((end-plotBackIndex)/2)),'b');
hold on; plot(Xcc,V(:,end-plotBackIndex),'r'); grid on;
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
xlim([0 xp1]);
%
subplot(2,3,3);
hold on; plot(Xcc,P(:,1),'black'); box on;
hold on; plot(Xcc,P(:,round((end-plotBackIndex)/2)),'b');
hold on; plot(Xcc,P(:,end-plotBackIndex),'r'); grid on;
hold on; plot(Xcc,Ti(:,end-plotBackIndex),'g');
hold on; plot(Xcc,Te(:,end-plotBackIndex),'g--');
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('P');
title('thermal pressure'); axis('square');
xlim([0 xp1]);
%
%
subplot(2,3,4);
hold on; plot(Xce,Ez(:,1),'black'); box on;
hold on; plot(Xce,Ez(:,round((end-plotBackIndex)/2)),'b');
hold on; plot(Xce,Ez(:,end-plotBackIndex),'r'); grid on;
hold on; plot(Xcc,Ez0(:,end-plotBackIndex),'g--'); 
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Ez');
title('electric field'); axis('square');
xlim([0 xp1]);
%
subplot(2,3,5);
hold on; plot(Xcc,B(:,1),'black'); box on;
hold on; plot(Xcc,B(:,round((end-plotBackIndex)/2)),'b'); grid on;
hold on; plot(Xcc,B(:,end-plotBackIndex),'r');
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field');
%set(gca,'xtick',0.0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 xp1]);
%
subplot(2,3,6);
hold on; plot(Xce,J(:,1),'black'); box on;
hold on; plot(Xce,J(:,round((end-plotBackIndex)/2)),'b'); grid on;
hold on; plot(Xce,J(:,end-plotBackIndex),'r');
%hold on; plot(Xce,J0(:,end-plotBackIndex),'g--');
hold on; plot(Xcc,Jcc(:,end-plotBackIndex),'g--');
xlabel('x'); ylabel('J'); axis('square');
title('current density');
%set(gca,'xtick',0:0.05:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 xp1]);
end



% figure(8); plot(tout,Poynt);
% xlabel('time'); title('poynting flux');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%    plot consevation laws
%%%
%%%


f8=figure(8); set(f8,'position',[1000 560 1400 800]);


%%%   plot mass conservation
%
subplot(2,3,1);
hold on; plot(tout,100*abs(Mass-Mass(1))/Mass(1)); box on;
set(gca,'yscale','log');
xlabel('time'); ylabel('% error'); title('mass conservation');
lg9 = legend('show'); set(lg9,'location','best');
xlim([0 tout(end)]);


%%%   plot momentum conservation
%
subplot(2,3,2);
MomGain = cumtrapz(tout,Msource);
%MomGain2 = cumtrapz(tout,Mflux);
hold on; plot(tout,Momentum,'displayName','\int NU dX'); box on;
hold on; plot(tout,MomGain,'linestyle','--','displayName','time-integrated source');
xlabel('time'); ylabel('momentum'); title('momentum conservation');
lg9 = legend('show'); set(lg9,'location','best');
xlim([0 tout(end)]);


%%%   plot entropy conservation 
%%%   (note most entropy comes from viscosity in shock front )
%
subplot(2,3,3);
EntropyGain = cumtrapz(tout,Ssource);
hold on; plot(tout,Entropy-Entropy(1),'displayName','\int S dX'); box on;
hold on; plot(tout,EntropyGain,'linestyle','--','displayName','time-integrated \eta J^2');
xlabel('time'); ylabel('entropy'); title('entropy conservation');
lg9 = legend('show'); set(lg9,'location','best');
xlim([0 tout(end)]);


%%%   plot energy conservation
%
subplot(2,3,5);
Etot2 = cumtrapz(tout,Poynt);
hold on; plot(tout,Etot-Etot(1),'displayName','energy'); box on;
hold on; plot(tout,Etot2,'linestyle','--','displayName','\int poynting flux');
xlabel('time'); ylabel('energy'); title('energy conservation');
lg9 = legend('show'); set(lg9,'location','northwest');
xlim([0 tout(end)]);


%%%   plot total magnetic field conservation
%
Btot2 = cumtrapz(tout,Bflux);
subplot(2,3,4); 
hold on; plot(tout,Btot,'displayName','\int B dx'); box on;
hold on; plot(tout,Btot2,'linestyle','--','displayName','time-integraged flux');
xlabel('time'); ylabel('magnetic flux'); title('magnetic flux conservation');
lg9 = legend('show'); set(lg9,'location','northwest');
xlim([0 tout(end)]);


subplot(2,3,6);
EthermGain = cumtrapz(tout,Psource);
hold on; plot(tout,Etherm-Etherm(1),'displayName','thermal'); box on;
hold on; plot(tout,Emean,'displayName','mean');
hold on; plot(tout,Efield,'displayName','magnetic'); box on;

xlabel('time'); ylabel('energy'); title('energy partition');
lg9 = legend('show'); set(lg9,'location','northwest');
xlim([0 tout(end)]);



%%%   plot JdotE and energy in gas
%
Jheating = cumtrapz(tout,JdotE);
f12= figure(12); 
hold on; plot(tout,Emean+Etherm-Etherm(1),'displayName', '\int (3/2P + \rhoV^2/2) dx');
hold on; plot(tout,Jheating,'linestyle','--','displayName','\int \int J\cdot E dx dt');
xlabel('time'); ylabel('energy'); title('kinetic energy conservation');
lg9 = legend('show'); set(lg9,'location','northwest');
box on; grid on;
xlim([0 tout(end)]);



%%%   calculate current-channel and shock location/speeds
%
X_J = zeros(size(tout));
X_N = zeros(size(tout));
X_shock = zeros(size(tout));
for n=1:length(tout)
    [~,ix0] = max(abs(J(:,n)));
    X_J(n) = Xce(ix0);
    %
    [~,ix0] = max(abs(N(:,n)));
    X_N(n) = Xcc(ix0);
    thisV = V(ix0,n);
    %
    if(thisV>0)
        [~,ix1] = min(abs(N(ix0:end,n)-2.0));
        ix1 = ix1 + ix0 - 1;
    else
        [~,ix1] = min(abs(N(1:ix0,n)-2.0));
        %ix1 = ix1 + ix0 - 1;
    end
    X_shock(n) = Xcc(ix1);
end

V_J = zeros(size(tout));
V_N = zeros(size(tout));
V_shock = zeros(size(tout));

for n=2:length(tout)
    V_J(n) = (X_J(n)-X_J(n-1))/(tout(n)-tout(n-1));
    V_N(n) = (X_N(n)-X_N(n-1))/(tout(n)-tout(n-1));
    V_shock(n) = (X_shock(n)-X_shock(n-1))/(tout(n)-tout(n-1));
end

[~,it0] = min(abs(tout-0.02));
[~,it1] = min(abs(tout-0.04));
V_J0 = (X_J(it1)-X_J(it0))/(0.04-0.02);
V_N0 = (X_N(it1)-X_N(it0))/(0.04-0.02);
V_shock0 = (X_shock(it1)-X_shock(it0))/(0.04-0.02);

f13=figure(13); 
plot(tout,X_J,tout,X_N,tout,X_shock); grid on;
xlabel('time'); ylabel('position');
title('axial positions');
lg13=legend('peak current density','peak density','shock front');
set(lg13,'location','best');
xlim([0 tout(end)]);



display(TeTi0);


