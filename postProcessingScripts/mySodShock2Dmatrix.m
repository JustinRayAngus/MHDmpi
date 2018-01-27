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

numProcs = 1;
filePath = '../physicsMods/sodShock2Dmatrix/';
%filePath = '../physicsMods/testing/';

fluxDir = 1; % 0 for x and 1 for z

for i=1:numProcs
fileName = ['output',num2str(i-1),'.h5'];
thisFile = [filePath,fileName];
procID  = hdf5read(thisFile,'procID');
fileinfo = hdf5info(thisFile);
Xcc = hdf5read(thisFile,'Xcc');
Xce = hdf5read(thisFile,'Xce');
Zcc = hdf5read(thisFile,'Zcc');
Zce = hdf5read(thisFile,'Zce');

N  = hdf5read(thisFile,'N');
Mx = hdf5read(thisFile,'Mx');
Mz = hdf5read(thisFile,'Mz');
E  = hdf5read(thisFile,'E');
P  = hdf5read(thisFile,'P');
Vx = hdf5read(thisFile,'Vx');
Vz = hdf5read(thisFile,'Vz');
Cs  = hdf5read(thisFile,'Cs');
gamma0 = hdf5read(thisFile,'gamma0');
%FluxLimLx = hdf5read(thisFile,'FluxLimLx');
%FluxLimRx    = hdf5read(thisFile,'FluxLimRx');
%FluxLx    = hdf5read(thisFile,'FluxLx');
%FluxRx    = hdf5read(thisFile,'FluxRx');
%FluxN_x  = hdf5read(thisFile,'FluxN_x');
FluxMx_x = hdf5read(thisFile,'FluxMx_x');
FluxMx_z = hdf5read(thisFile,'FluxMx_z');
FluxMz_x = hdf5read(thisFile,'FluxMz_x');
FluxMz_z = hdf5read(thisFile,'FluxMz_z');
%FluxMz_x = hdf5read(thisFile,'FluxMz_x');
%FluxE_x  = hdf5read(thisFile,'FluxE_x');
tout= hdf5read(thisFile,'tout');



f11=figure(11); 
%set(f1,'position',[1030 925 1100 420]);
set(f11,'position',[1800 360 500 760]);
subplot(2,1,1);
hold on; surf(Xcc,Zcc,N(:,:,1)); shading flat;colorbar;
xlabel('x direction');  %caxis([0,1]);
ylabel('z direction');
title('initial density profile');
axis([-0.5 0.5 -0.5 0.5]); axis('equal');
%
subplot(2,1,2);
hold on; surf(Xcc,Zcc,N(:,:,round(end))); colorbar;
xlabel('x direction');  %caxis([0,1]);
ylabel('z direction'); shading flat;
title('density profile at final time step');
axis([-0.5 0.5 -0.5 0.5]); axis('equal');
%
map = colormap('jet');
map(1,:) = [1 1 1];
colormap(map)



%%%
%

f1=figure(1); 
set(f1,'position',[540 1 900 800]);
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



end

