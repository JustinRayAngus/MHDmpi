%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 2D pinch module of m=0 mode
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


colormap('jet');
%colormap('hot');

filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/testing7/'; numProcs = 20; newDeck = 2;





plotBackIndex = 1; % plot time will be end-plotBackIndex
thist = 10; %9.4; %7.9;
thist2 = 14.2; %10.2; %9.4;

tout = loadData(filePath,numProcs,'tout');
%tout = tout*t0/tA; % normalize to 2010 paper alfven time
[~,tindex] = min(abs(tout-thist2));
[~,tindex2] = min(abs(tout-tout(end-1)));

Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Zcc = loadData(filePath,numProcs,'Zcc');
Zce = loadData(filePath,numProcs,'Zce');
rcc = loadData(filePath,numProcs,'rcc');

N  = loadData(filePath,numProcs,'N');
deltaP  = loadData(filePath,numProcs,'deltaP');
%P0  = hdf5read(thisFile,'P0');
Mx  = loadData(filePath,numProcs,'Mx');
Mz  = loadData(filePath,numProcs,'Mz');
S  = loadData(filePath,numProcs,'S');
By = loadData(filePath,numProcs,'By');
P  = loadData(filePath,numProcs,'P');
P0  = loadData(filePath,numProcs,'P0');

Vx  = loadData(filePath,numProcs,'Vx');
Vz  = loadData(filePath,numProcs,'Vz');
Jz  = loadData(filePath,numProcs,'Jz');
Jx  = loadData(filePath,numProcs,'Jx');
Ez  = loadData(filePath,numProcs,'Ez');
Ex  = loadData(filePath,numProcs,'Ex');
Cs  = loadData(filePath,numProcs,'Cs');
gamma0 = loadData(filePath,numProcs,'gamma0');
dX = Xcc(2)-Xcc(1);
dZ = Zcc(2)-Zcc(1);
Fx = loadData(filePath,numProcs,'Fx'); % - d(P+B^2/2)/dr - B^2/r

Jz0 = Jz;
Ez0 = Ez;
if(newDeck>=1)
  Jz0  = loadData(filePath,numProcs,'Jz0');
  Ez0  = loadData(filePath,numProcs,'Ez0');
  if(newDeck>=2)
      delta0 = loadData(filePath,numProcs,'delta0');
      lambda0 = loadData(filePath,numProcs,'lambda0');
   %   Ee_source = loadData(filePath,numProcs,'Ee_source');
   %   Ei_source = loadData(filePath,numProcs,'Ei_source');
      Te = loadData(filePath,numProcs,'Te');
      Ti = loadData(filePath,numProcs,'Ti');
      Vex = Vx - lambda0*Jx./N;
      Vez = Vz - lambda0*Jz./N;
   %   Ee  = loadData(filePath,numProcs,'Ee');
   %   Ei  = loadData(filePath,numProcs,'Ei');
  end
end


T = P./N;
Ptot = P+By.^2/2+Mx.*Vx/2.0;
%figure(2); hold on; plot(Xce,FluxR(:,5),'r--');
%figure(4); hold on; plot(Xce,FluxLim(:,end),'b');
%figure(4); hold on; plot(Xce,FluxRatio(:,end),'r');
%%%

fluxDir = 0;
if(fluxDir==0)
   index0 = round(length(Zcc)/2); 
   index1 = 1:length(Xcc);
   xlab = 'x';
   plotVec = Xcc;
   VelVec = Vx;
end



plots=1;
if(plots)
f1=figure(7); 
set(f1,'position',[1030 425 1300 840]);
%set(f1,'position',[341 436 900 840]);

subplot(2,3,1);
hold on; plot(Xcc,Fx(index0,:,1),'black'); box on;
hold on; plot(Xcc,Fx(index0,:,2),'r'); grid on;
hold on; plot(Xcc,Fx(index0,:,20),'g');
hold on; plot(Xcc,Fx(index0,:,80),'cyan');
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Fx');
title('radial force'); axis('square');
xlim([0 Xce(end-1)]);
%
subplot(2,3,2);
hold on; plot(Xcc,Vx(index0,index1,1),'black'); box on;
hold on; plot(Xcc,Vx(index0,index1,2),'r'); grid on;
hold on; plot(Xcc,Vx(index0,:,20),'g');
hold on; plot(Xcc,Vx(index0,:,80),'cyan');
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
xlim([0 Xce(end-1)]);
%
subplot(2,3,3);
hold on; plot(Xcc,Ez(index0,:,1),'black'); box on;
hold on; plot(Xcc,Ez(index0,:,2),'r'); grid on;
hold on; plot(Xcc,Ez(index0,:,20),'g'); grid on;
hold on; plot(Xcc,Ez(index0,:,60),'b'); grid on;
hold on; plot(Xcc,Ez(index0,:,80),'cyan'); grid on;
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Ez');
title('electric field'); axis('square');
xlim([0 Xce(end-1)]);
%
%
subplot(2,3,4);
hold on; plot(Xcc,P(index0,index1,1),'black'); box on;
hold on; plot(Xcc,P(index0,index1,2),'r'); grid on;
hold on; plot(Xcc,P(index0,index1,20),'g');
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('P');
title('pressure'); axis('square');
xlim([0 Xce(end-1)]);
ylim([0 1.2]);
%
subplot(2,3,5);
hold on; plot(Xcc,By(index0,index1,1),'black'); box on;
hold on; plot(Xcc,By(index0,index1,2),'r');
hold on; plot(Xcc,By(index0,index1,20),'g');
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field'); grid on;
%set(gca,'xtick',0.0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 Xce(end-1)]);
%
subplot(2,3,6);
hold on; plot(Xcc,Jz(index0,:,1),'black'); box on; grid on;
hold on; plot(Xcc,Jz(index0,:,2),'r');
hold on; plot(Xcc,Jz(index0,:,20),'g');
xlabel('x'); ylabel('J_z'); axis('square');
title('current density');
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 Xce(end-1)]);
end






