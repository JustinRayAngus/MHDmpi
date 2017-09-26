%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 1D Sod shock simulations
%%%
%%%   see G.A. Sod,. J. Comp. Phys. 27, 1-31 (1978)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

numProcs = 4;
filePath = '../physicsMods/sodShock/';

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
%close(figure(1));
f1=figure(1); 
set(f1,'position',[1030 925 1100 420]);
%set(f1,'position',[341 436 1100 420]);


subplot(1,2,1);
hold on; plot(Xcc,N(:,1)); box on;
hold on; plot(Xcc,N(:,round(end/2)));
hold on; plot(Xcc,N(:,end)); grid on;
set(gca,'xtick',-0.5:0.1:0.5);
xlabel('x'); ylabel('f');
title('mass density');
xlim([-0.5 0.5]);
%
subplot(1,2,2);
hold on; plot(Xce,FluxN(:,1)); box on;
hold on; plot(Xce,FluxN(:,round(end/2)));
hold on; plot(Xce,FluxN(:,round(end))); grid on;
set(gca,'xtick',-0.5:0.1:0.5);
xlabel('x'); ylabel('flux');
title('flux');
xlim([-0.5 0.5]);

%
%
%

f2=figure(2); 
set(f2,'position',[1030 425 1100 420]);
%set(f1,'position',[341 436 1100 420]);


subplot(1,2,1);
hold on; plot(Xcc,P(:,1)); box on;
hold on; plot(Xcc,P(:,round(end/2)));
hold on; plot(Xcc,P(:,end)); grid on;
set(gca,'xtick',-0.5:0.1:0.5);
xlabel('x'); ylabel('f');
title('momentum density');
xlim([-0.5 0.5]);
%
subplot(1,2,2);
hold on; plot(Xce,FluxM(:,1)); box on;
hold on; plot(Xce,FluxM(:,round(end/2)));
hold on; plot(Xce,FluxM(:,round(end))); grid on;
set(gca,'xtick',-0.5:0.1:0.5);
xlabel('x'); ylabel('flux');
title('flux');
xlim([-0.5 0.5]);


figure(4);
hold on; plot(Xcc,P(:,end),'black'); box on;
hold on; plot(Xcc,N(:,end),'black--'); box on;
axis([-0.5 0.5 0 1.2]); grid on;
xlabel('x'); ylabel('density (solid), pressure (dashed)');
title(['t=',num2str(tout(end)),': compare with sodShock BOUT++']);

end

