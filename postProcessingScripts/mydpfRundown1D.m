%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 1D Sod shock simulations
%%%
%%%   see G.A. Sod,. J. Comp. Phys. 27, 1-31 (1978)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

numProcs = 4;
filePath = '../physicsMods/dpfRundown1D/';


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
J0  = hdf5read(thisFile,'J0'); % curl(B)
Ez  = hdf5read(thisFile,'Ez');
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
FluxEz  = hdf5read(thisFile,'FluxEz');
tout= hdf5read(thisFile,'tout');
dX = Xcc(2)-Xcc(1);


T = P./N;
eta0 = 1e-3;
figure(2); hold on; plot(Xcc,P(:,2),'black*');
%P2 = S.*N.^(5/3-1);
%figure(2); hold on; plot(Xce,FluxR(:,5),'r--');
figure(4); hold on; plot(Xce,FluxLim(:,2),'b');
figure(4); hold on; plot(Xce,FluxRatio(:,2),'r');
%%%
%

plots=1;
if(plots)
f1=figure(1); 
set(f1,'position',[1030 425 1300 840]);
%set(f1,'position',[341 436 900 840]);

subplot(2,3,1);
hold on; plot(Xcc,N(:,1),'black'); box on;
hold on; plot(Xcc,N(:,round(end/2)),'b');
hold on; plot(Xcc,N(:,end),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('N');
title('mass density'); axis('square');
xlim([-0.51 0.5]);
%
subplot(2,3,2);
hold on; plot(Xcc,M(:,1),'black'); box on;
hold on; plot(Xcc,M(:,round(end/2)),'b');
hold on; plot(Xcc,M(:,end),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
xlim([-0.51 0.5]);
%
subplot(2,3,3);
hold on; plot(Xcc,P(:,1),'black'); box on;
hold on; plot(Xcc,P(:,round(end/2)),'b');
hold on; plot(Xcc,P(:,end),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('P');
title('thermal pressure'); axis('square');
xlim([-0.51 0.5]);
%
%
subplot(2,3,4);
hold on; plot(Xcc,Ez(:,1),'black'); box on;
hold on; plot(Xcc,Ez(:,round(end/2)),'b');
hold on; plot(Xcc,Ez(:,end),'r'); grid on;
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Ez');
title('electric field'); axis('square');
xlim([-0.51 0.5]);
%
subplot(2,3,5);
eta = P./N/(gamma0-1);
hold on; plot(Xcc,B(:,1),'black'); box on;
hold on; plot(Xcc,B(:,round(end/2)),'b'); grid on;
hold on; plot(Xcc,B(:,end),'r');
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field');
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',1:0.5:3);
xlim([-0.5 0.5]);
%
subplot(2,3,6);
hold on; plot(Xcc,J(:,1),'black'); box on;
hold on; plot(Xcc,J(:,round(end/2)),'b'); grid on;
hold on; plot(Xcc,J(:,end),'r');
hold on; plot(Xcc,J0(:,end),'g--');
xlabel('x'); ylabel('J'); axis('square');
title('current density');
set(gca,'xtick',-0.5:0.25:0.5);
%set(gca,'ytick',1:0.5:3);
axis([-0.5 0.5 -10 0]);
end


end

