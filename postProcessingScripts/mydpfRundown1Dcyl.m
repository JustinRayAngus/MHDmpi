%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D dpf rundown in cyl coords
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
filePath = '../physicsMods/dpfRundown1Dcyl/';

plotBackIndex = 0; % plot time will be end-plotBackIndex
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
S  = loadData1DVec(filePath,numProcs,'S');
B  = loadData1DVec(filePath,numProcs,'B');
P  = loadData1DVec(filePath,numProcs,'P');
T  = loadData1DVec(filePath,numProcs,'T');
V  = loadData1DVec(filePath,numProcs,'V');
J  = loadData1DVec(filePath,numProcs,'J');
Jcc = loadData1DVec(filePath,numProcs,'Jcc');
J0  = loadData1DVec(filePath,numProcs,'J0');
Ez  = loadData1DVec(filePath,numProcs,'Ez');
eta = loadData1DVec(filePath,numProcs,'eta');
etace = loadData1DVec(filePath,numProcs,'etace');
gamma0 = loadData1DVec(filePath,numProcs,'gamma0');
delta0 = loadData1DVec(filePath,numProcs,'delta0');
%
FluxEz  = loadData1DVec(filePath,numProcs,'FluxEz');
FluxB  = loadData1DVec(filePath,numProcs,'FluxB');
FluxN  = loadData1DVec(filePath,numProcs,'FluxN');
FluxS  = loadData1DVec(filePath,numProcs,'FluxS');
%
rcc  = loadData1DVec(filePath,numProcs,'rcc');
rce  = loadData1DVec(filePath,numProcs,'rce');

%%%   calculate divU
%
dX = Xcc(2)-Xcc(1);
divV = zeros(size(V));
dPdx = zeros(size(V));
for j=2:length(Xcc)-1
    divV(j,:) = (V(j+1,:)-V(j-1,:))/2/dX;
    dPdx(j,:) = (P(j+1,:)-P(j-1,:))/2/dX;
end


%%%   calcualte Ez at cell-center
%
Ezcc = zeros(size(N));
for j=3:length(Xcc)-2
    Ezcc(j,:) = (Ez(j-1,:) + Ez(j,:))/2;
end


%figure(77); hold on; plot(Xcc,divV(:,100)); box on;

Ptot = P + B.^2/2 + M.*V/2.0;
PdV = P.*divV;
VdP = V.*dPdx;
etaJ2 = eta.*Jcc.^2;
etaJ2mod = (gamma0-1)*etaJ2./N.^(gamma0-1);
Edenstot = 3/2*P + N.*V.^2/2 + B.^2/2;

%%%   calculate conservation stuff
%
dV     = 2*pi*Xcc(3:end-2)*dX;
JdotE  = sum(Jcc(3:end-2,:).*Ezcc(3:end-2,:))*dX;
Etot   = sum(Edenstot(3:end-2,:))*dX;
Etherm = sum(1.5*P(3:end-2,:))*dX;
Emean  = sum(0.5*M(3:end-2,:).*V(3:end-2,:))*dX;
Efield = sum(B(3:end-2,:).^2/2)*dX;
EEztot = delta0*sum(Ez(2:end-2,:).^2/2.0)*dX;
%
Entropy  = sum(S(3:end-2,:).*dV);
Momentum = sum(M(3:end-2,:).*dV);
Mass     = sum(N(3:end-2,:).*dV);

for n=1:length(tout)

    Btot(n)   = sum(B(3:end-2,n))*dX;
    Ssource(n) = sum(etaJ2mod(3:end-2,n))*dX;
    Msource(n) = sum(-Jcc(3:end-2,n).*B(3:end-2,n))*dX;
    Psource(n) = sum(etaJ2(3:end-2,n)-0*PdV(3:end-2,n))*dX;
   % Psource(n) = sum(etaJ2(3:end-2,n)+VdP(3:end-2,n))*dX;
    Etot(n)    = Etot(n) + EEztot(n);
    
    Poynt(n) = -Ez(2,n)*B(3,n) ...
             - 0*5/2*(P(2,n)+P(3,n))/2*(V(2,n)+V(3,n))/2 ...
             - 0*0.5*N(3,n)*V(3,n)^3;
    Bflux(n) = -Ez(2,n);
    Mflux(n) =  (P(2,n)+P(3,n))/2 + B(2,n)^2/2;
    Msource(n) = Msource(n) + (P(2,n)+P(3,n))/2;
  %  Psource(n) = Psource(n) - 3/2*P(2,n)*(V(2,n)+V(3,n))/2;

    Poynt(n) = Poynt(n) + Ez(end-1,n)*B(end-2,n) + 5/2*P(end-1,n)*V(end-1,n);
    Msource(n) = Msource(n) - (P(end-2,n)+P(end-1,n))/2;
    Mflux(n) = Mflux(n) - (P(end-2,n)+P(end-1,n))/2;
   
end
Cs = sqrt(gamma0*P./N);
Mach = abs(V./Cs);
Cspeed = abs(V) + sqrt(gamma0*P./N + B.^2./N);


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
hold on; plot(Xcc,N(:,round(end/2)),'b');
%hold on; plot(Xcc,N(:,round(120)),'g');
%hold on; plot(Xcc,N(:,round(end/2)),'b');
hold on; plot(Xcc,N(:,end-plotBackIndex),'r'); grid on;
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('r'); ylabel('N');
title('mass density'); axis('square');
xlim([0 xp1]);
%
subplot(2,3,2);
hold on; plot(Xcc,V(:,1),'black'); box on;
hold on; plot(Xcc,V(:,round(end/2)),'b');
hold on; plot(Xcc,V(:,end-plotBackIndex),'r'); grid on;
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('r'); ylabel('V');
title('velocity'); axis('square');
xlim([0 xp1]);
%
subplot(2,3,3);
hold on; plot(Xcc,P(:,1),'black'); box on;
hold on; plot(Xcc,P(:,round(end/2)),'b');
hold on; plot(Xcc,P(:,end-plotBackIndex),'r'); grid on;
hold on; plot(Xcc,T(:,end-plotBackIndex),'g');
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('r'); ylabel('P');
title('thermal pressure'); axis('square');
xlim([0 xp1]);
%
%
subplot(2,3,4);
hold on; plot(Xce,Ez(:,1),'black'); box on;
hold on; plot(Xce,Ez(:,round(end/2)),'b');
hold on; plot(Xce,Ez(:,end-plotBackIndex),'r'); grid on;
hold on; plot(Xcc,Ez0(:,end-plotBackIndex),'g--'); 
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('r'); ylabel('Ez');
title('electric field'); axis('square');
xlim([0 xp1]);
%
subplot(2,3,5);
hold on; plot(Xcc,B(:,1),'black'); box on;
hold on; plot(Xcc,B(:,round(end/2)),'b'); grid on;
hold on; plot(Xcc,B(:,end-plotBackIndex),'r');
xlabel('r'); ylabel('B'); axis('square');
title('magnetic field');
%set(gca,'xtick',0.0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 xp1]);
%
subplot(2,3,6);
hold on; plot(Xce,J(:,1),'black'); box on;
hold on; plot(Xce,J(:,round(end/2)),'b'); grid on;
hold on; plot(Xce,J(:,end-plotBackIndex),'r');
hold on; plot(Xce,J0(:,end-plotBackIndex),'g--');
%hold on; plot(Xcc,Jcc(:,end-plotBackIndex),'g--');
xlabel('r'); ylabel('J'); axis('square');
title('current density');
%set(gca,'xtick',0:0.05:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 xp1]);
end

figure(7);
hold on; plot(Xcc,Mach(:,end-plotBackIndex),'r'); box on;
title('Mach Number'); xlim([0 xp1]);



% figure(8); plot(tout,Poynt);
% xlabel('time'); title('poynting flux');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%    plot consevation laws
%%%
%%%


% f8=figure(8); set(f8,'position',[1000 560 1400 800]);


% %%%   plot mass conservation
% %
% subplot(2,3,1);
% hold on; plot(tout,(Mass-Mass(1))/Mass(1)); box on;
% xlabel('time'); ylabel('% error'); title('mass conservation');
% lg9 = legend('show'); set(lg9,'location','best');
% xlim([0 tout(end)]);
% 
% 
% %%%   plot momentum conservation
% %
% subplot(2,3,2);
% MomGain = cumtrapz(tout,Msource);
% %MomGain2 = cumtrapz(tout,Mflux);
% hold on; plot(tout,Momentum,'displayName','\int NU dX'); box on;
% hold on; plot(tout,MomGain,'linestyle','--','displayName','time-integrated source');
% xlabel('time'); ylabel('momentum'); title('momentum conservation');
% lg9 = legend('show'); set(lg9,'location','best');
% xlim([0 tout(end)]);
% 
% 
% %%%   plot entropy conservation
% %
% subplot(2,3,3);
% EntropyGain = cumtrapz(tout,Ssource);
% hold on; plot(tout,Entropy-Entropy(1),'displayName','\int S dX'); box on;
% hold on; plot(tout,EntropyGain,'linestyle','--','displayName','time-integrated source');
% xlabel('time'); ylabel('entropy'); title('entropy conservation');
% lg9 = legend('show'); set(lg9,'location','best');
% xlim([0 tout(end)]);
% 
% 
% %%%   plot energy conservation
% %
% subplot(2,3,5);
% Etot2 = cumtrapz(tout,Poynt);
% hold on; plot(tout,Etot-Etot(1),'displayName','energy'); box on;
% hold on; plot(tout,Etot2,'linestyle','--','displayName','\int poynting flux');
% xlabel('time'); ylabel('energy'); title('energy conservation');
% lg9 = legend('show'); set(lg9,'location','northwest');
% xlim([0 tout(end)]);
% 
% 
% %%%   plot total magnetic field conservation
% %
% Btot2 = cumtrapz(tout,Bflux);
% subplot(2,3,4); 
% hold on; plot(tout,Btot,'displayName','\int B dx'); box on;
% hold on; plot(tout,Btot2,'linestyle','--','displayName','time-integraged flux');
% xlabel('time'); ylabel('magnetic flux'); title('magnetic flux conservation');
% lg9 = legend('show'); set(lg9,'location','northwest');
% xlim([0 tout(end)]);
% 
% 
% subplot(2,3,6);
% EthermGain = cumtrapz(tout,Psource);
% hold on; plot(tout,Etherm-Etherm(1),'displayName','thermal energy'); box on;
% hold on; plot(tout,EthermGain,'linestyle','--','displayName','\int \eta J^2-P\nabla U flux');
% xlabel('time'); ylabel('energy'); title('thermal energy conservation');
% lg9 = legend('show'); set(lg9,'location','northwest');
% xlim([0 tout(end)]);
% 
% 
% 
% %%%   plot JdotE and energy in gas
% %
% Jheating = cumtrapz(tout,JdotE);
% figure(12); 
% hold on; plot(tout,Emean+Etherm-Etherm(1),'displayName', '\int (3/2P + \rhoV^2/2) dx');
% hold on; plot(tout,Jheating,'displayName','\int \int J\cdot E dx dt');
% xlabel('time'); ylabel('energy'); title('kinetic energy conservation');
% lg9 = legend('show'); set(lg9,'location','northwest');
% box on; grid on;
% xlim([0 tout(end)]);

