%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 1D Sod shock simulations
%%%
%%%   see G.A. Sod,. J. Comp. Phys. 27, 1-31 (1978)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

numProcs = 4;
filePath = '../physicsMods/mhdShock/';


for i=1:numProcs
fileName = ['output',num2str(i-1),'.h5'];
thisFile = [filePath,fileName];
procID  = hdf5read(thisFile,'procID');
fileinfo = hdf5info(thisFile);
Xcc = hdf5read(thisFile,'Xcc');
Xce = hdf5read(thisFile,'Xce');
N  = hdf5read(thisFile,'N');
M  = hdf5read(thisFile,'M');
S  = hdf5read(thisFile,'S');
B  = hdf5read(thisFile,'B');
P  = hdf5read(thisFile,'P');
V  = hdf5read(thisFile,'V');
J  = hdf5read(thisFile,'J');
Cs  = hdf5read(thisFile,'Cs');
gamma0 = hdf5read(thisFile,'gamma0');
FluxRatio  = hdf5read(thisFile,'FluxRatio');
FluxLim    = hdf5read(thisFile,'FluxLim');
FluxL    = hdf5read(thisFile,'FluxL');
FluxR    = hdf5read(thisFile,'FluxR');
FluxN  = hdf5read(thisFile,'FluxN');
FluxM  = hdf5read(thisFile,'FluxM');
FluxS  = hdf5read(thisFile,'FluxS');
FluxB  = hdf5read(thisFile,'FluxB');
tout= hdf5read(thisFile,'tout');
dX = Xcc(2)-Xcc(1);

figure(2); hold on; plot(Xcc,J(:,1));

%%%
%

f1=figure(1); 
set(f1,'position',[1030 425 900 840]);
%set(f1,'position',[341 436 900 840]);

subplot(2,2,1);
hold on; plot(Xcc,N(:,1),'black'); box on;
hold on; plot(Xcc,N(:,4),'b');
hold on; plot(Xcc,N(:,8),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('N');
title('mass density'); axis('square');
xlim([-0.51 0.5]);
%
subplot(2,2,2);
hold on; plot(Xcc,V(:,1),'black'); box on;
hold on; plot(Xcc,V(:,4),'b');
hold on; plot(Xcc,V(:,8),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
xlim([-0.51 0.5]);
%
subplot(2,2,3);
hold on; plot(Xcc,P(:,1),'black'); box on;
hold on; plot(Xcc,P(:,4),'b');
hold on; plot(Xcc,P(:,6),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('P');
title('thermal pressure'); axis('square');
xlim([-0.51 0.5]);
%
subplot(2,2,4);
eta = P./N/(gamma0-1);
hold on; plot(Xcc,B(:,1),'black'); box on;
hold on; plot(Xcc,B(:,4),'b'); grid on;
hold on; plot(Xcc,B(:,8),'r');
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field');
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',1:0.5:3);
xlim([-0.51 0.5]);


end

