%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 1D Sod shock simulations
%%%
%%%   see G.A. Sod,. J. Comp. Phys. 27, 1-31 (1978)
%%%
%%%   Van Leer works better than suber bee 
%%%   see "A Primer on Eulerian Computational Fluid Dynamics for Astrophysics"
%%%   by Hy Trac and Ue-Li Pen pg 309 (2003)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'defaultaxeslinewidth',1.5);
set(0,'defaultlinelinewidth',1.5);
set(0,'defaultaxesfontweight','bold');


numProcs = 4;
filePath = '../physicsMods/sodShock/data_nx200/';
%filePath = '../physicsMods/sodShock/data_nx400/';
filePath = '../physicsMods/sodShock/data_nx800/';
filePath = '../physicsMods/sodShock/data_nx800/';

for i=1:numProcs
fileName = ['output',num2str(i-1),'.h5'];
thisFile = [filePath,fileName];
procID  = hdf5read(thisFile,'procID');
fileinfo = hdf5info(thisFile);
Xcc = hdf5read(thisFile,'Xcc');
Xce = hdf5read(thisFile,'Xce');
N  = hdf5read(thisFile,'N');
M  = hdf5read(thisFile,'M');
E  = hdf5read(thisFile,'E');
P  = hdf5read(thisFile,'P');
V  = hdf5read(thisFile,'V');
Cs  = hdf5read(thisFile,'Cs');
gamma0 = hdf5read(thisFile,'gamma0');
FluxRatio  = hdf5read(thisFile,'FluxRatio');
FluxLim    = hdf5read(thisFile,'FluxLim');
FluxL    = hdf5read(thisFile,'FluxL');
FluxR    = hdf5read(thisFile,'FluxR');
FluxN  = hdf5read(thisFile,'FluxN');
FluxM  = hdf5read(thisFile,'FluxM');
FluxE  = hdf5read(thisFile,'FluxE');
tout= hdf5read(thisFile,'tout');


%%%
%

f2=figure(22); 
set(f2,'position',[540 1 900 800]);
%set(f1,'position',[341 436 900 840]);

subplot(2,2,1);
hold on; plot(Xcc,N(:,1),'black'); box on;
hold on; plot(Xcc,N(:,4),'b');
hold on; plot(Xcc,N(:,8),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('N');
title('mass density'); axis('square');
axis([-0.5 0.5 0 1.2]);
%
subplot(2,2,2);
hold on; plot(Xcc,V(:,1),'black'); box on;
hold on; plot(Xcc,V(:,4),'b');
hold on; plot(Xcc,V(:,8),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
axis([-0.5 0.5 0 1.2]);
%
subplot(2,2,3);
hold on; plot(Xcc,P(:,1),'black'); box on;
hold on; plot(Xcc,P(:,4),'b');
hold on; plot(Xcc,P(:,8),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('P');
title('pressure'); axis('square');
axis([-0.5 0.5 0 1.2]);
%
subplot(2,2,4);
eta = P./N/(gamma0-1);
hold on; plot(Xcc,eta(:,1),'black'); box on;
hold on; plot(Xcc,eta(:,4),'b'); grid on;
hold on; plot(Xcc,eta(:,8),'r');
xlabel('x'); ylabel('\eta'); axis('square');
title('internal energy per unit mass');
set(gca,'xtick',-0.5:0.25:0.5);
set(gca,'ytick',1:0.5:3);
axis([-0.5 0.5 1 3]);


f4=figure(4); set(f4,'position', [1945 45 560 420]); 
hold on; plot(Xcc,P(:,end),'black'); box on;
hold on; plot(Xcc,N(:,end),'black--'); box on;
axis([-0.5 0.5 0 1.2]); grid on;
xlabel('x'); ylabel('pressure (solid), density (dashed)');
title(['t=',num2str(tout(end)),': compare with sodShock BOUT++']);


f5=figure(5); set(f5,'position', [1945 45 560 420]); 
hold on; plot(Xcc,P(:,2),'black'); box on;
hold on; plot(Xcc,P(:,3),'r'); box on;
axis([0 0.08 0.3 0.306]); grid on;
xlabel('x');
title('bad oscillations early on using U1 with nx>3000');


f6=figure(6); %set(f6,'position', [1945 45 560 420]); 
hold on; plot(Xce,FluxRatio(:,2),'black'); box on;
hold on; plot(Xce,FluxRatio(:,end),'r'); box on;
%axis([-0.5 0.5 0 1.2]); grid on;
xlabel('x');
title('flux ratio');
end

