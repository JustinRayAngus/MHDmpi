%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 2D Sod shock simulations using matrix2D
%%%
%%%   see G.A. Sod,. J. Comp. Phys. 27, 1-31 (1978)
%%%
%%%   Van Leer works better than suber bee 
%%%   see "A Primer on Eulerian Computational Fluid Dynamics for Astrophysics"
%%%   by Hy Trac and Ue-Li Pen pg 309 (2003)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

numProcs = 4;
filePath = '../physicsMods/sodShock2Dmatrix/';
filePath = '../physicsMods/sodShock2Dmatrix/data_weno5_refined/';

fluxDir = 0; % 0 for x and 1 for z


tout = loadData(filePath,numProcs,'tout');
Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Zcc = loadData(filePath,numProcs,'Zcc');
Zce = loadData(filePath,numProcs,'Zce');
%
N  = loadData(filePath,numProcs,'N');
Vx = loadData(filePath,numProcs,'Vx');
Vz = loadData(filePath,numProcs,'Vz');
P  = loadData(filePath,numProcs,'P');
E  = loadData(filePath,numProcs,'E');
Cs = loadData(filePath,numProcs,'Cs');
gamma0 = loadData(filePath,numProcs,'gamma0');



%fileName = ['output',num2str(i-1),'.h5'];
%thisFile = [filePath,fileName];
%procID  = hdf5read(thisFile,'procID');
%fileinfo = hdf5info(thisFile);


f11=figure(11); 
set(f11,'position',[1100 400 900 800]);
%set(f1,'position',[341 436 900 840]);

if(fluxDir==0)
   index0 = round(length(Zcc)/2); 
   index1 = 1:length(Xcc);
   xlab = 'x';
   plotVec = Xcc;
   VelVec = Vx;
end
if(fluxDir==1)
   index0 = 1:length(Zcc);
   index1 = round(length(Xcc)/2); 
   xlab = 'z';  
   plotVec = Zcc;
   VelVec = Vz;
end
subplot(2,2,1);
hold on; plot(plotVec,N(index0,index1,1),'black'); box on;
hold on; plot(plotVec,N(index0,index1,4),'b');
hold on; plot(plotVec,N(index0,index1,8),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
set(gca,'ytick',0:0.3:1.2);
xlabel(xlab); ylabel('N');
title('mass density'); axis('square');
axis([-0.5 0.5 0 1.2]);
%
subplot(2,2,2);
hold on; plot(plotVec,VelVec(index0,index1,1),'black'); box on;
hold on; plot(plotVec,VelVec(index0,index1,4),'b');
hold on; plot(plotVec,VelVec(index0,index1,8),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
set(gca,'ytick',0:0.3:1.2);
xlabel(xlab); ylabel('V');
title('velocity'); axis('square');
axis([-0.5 0.5 0 1.2]);
%
subplot(2,2,3);
hold on; plot(plotVec,P(index0,index1,1),'black'); box on;
hold on; plot(plotVec,P(index0,index1,4),'b');
hold on; plot(plotVec,P(index0,index1,8),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
set(gca,'ytick',0:0.3:1.2);
xlabel(xlab); ylabel('P');
title('pressure'); axis('square');
axis([-0.5 0.5 0 1.2]);
%
subplot(2,2,4);
eta = P./N/(gamma0-1);
hold on; plot(plotVec,eta(index0,index1,1),'black'); box on;
hold on; plot(plotVec,eta(index0,index1,4),'b'); grid on;
hold on; plot(plotVec,eta(index0,index1,8),'r');
xlabel(xlab); ylabel('\eta'); axis('square');
title('internal energy per unit mass');
set(gca,'xtick',-0.5:0.25:0.5);
set(gca,'ytick',1:0.5:3);
axis([-0.5 0.5 1 3]);


f4=figure(4); set(f4,'position', [1945 45 560 420]); 
hold on; plot(plotVec,P(index0,index1,end),'black'); box on;
hold on; plot(plotVec,N(index0,index1,end),'black--'); box on;
axis([-0.5 0.5 0 1.2]); grid on;
xlabel('x'); ylabel('pressure (solid), density (dashed)');
title(['t=',num2str(tout(end)),': compare with sodShock BOUT++']);


f5=figure(5); set(f5,'position', [1945 45 560 420]); 
hold on; plot(plotVec,P(index0,index1,2),'black'); box on;
hold on; plot(plotVec,P(index0,index1,3),'r'); box on;
axis([0 0.08 0.3 0.306]); grid on;
xlabel('x');
title('bad oscillations early on using U1 with nx>3000');


