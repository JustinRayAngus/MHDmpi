%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from MHD1D simulation
%%%
%%%   compare results with analytical solution from Burgers1DSolution.m
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

numProcs = 4;
filePath = '../physicsMods/burgers1D/';
%filePath = '../physicsMods/sodShock/';

for i=1:numProcs
fileName = ['output',num2str(i-1),'.h5'];
thisFile = [filePath,fileName];
procID  = hdf5read(thisFile,'procID');
fileinfo = hdf5info(thisFile);
Xcc = hdf5read(thisFile,'Xcc');
Xce = hdf5read(thisFile,'Xce');
F0  = hdf5read(thisFile,'F0');
FluxRatio  = hdf5read(thisFile,'FluxRatio');
FluxLim    = hdf5read(thisFile,'FluxLim');
FluxL    = hdf5read(thisFile,'FluxL');
FluxR    = hdf5read(thisFile,'FluxR');
Flux  = hdf5read(thisFile,'Flux');
tout= hdf5read(thisFile,'tout');


%%%
%
%close(figure(1));
f1=figure(1); 
%set(f1,'position',[1030 925 1100 420]);
set(f1,'position',[858 37 500 760]);


subplot(2,1,1);
hold on; plot(Xcc,F0(:,1)); box on;
hold on; plot(Xcc,F0(:,round(end/2)));
hold on; plot(Xcc,F0(:,end)); grid on;
set(gca,'xtick',-0.5:0.1:0.5);
xlabel('x'); ylabel('f');
title('function');
xlim([-0.5 0.5]);
%
subplot(2,1,2);
%hold on; plot(Xce,FluxL(:,1)); box on;
hold on; plot(Xce,FluxL(:,round(end/2)),'displayName','Left');
hold on; plot(Xce,FluxR(:,round(end/2)),'displayName','right');
%hold on; plot(Xce,FluxL(:,round(end))); grid on;
set(gca,'xtick',-0.5:0.1:0.5);
xlabel('x'); ylabel('flux'); grid on;
title('flux'); xlim([-0.5 0.5]);
if(i==1) 
    lg1=legend('show');
end

%
%
%

% f3=figure(3); set(f3,'position',[1030 430 1100 420]);
% %
% subplot(1,2,1);
% hold on; plot(Xce,FluxRatio(:,end)); grid on;
% set(gca,'xtick',-0.5:0.1:0.5); box on;
% xlim([-0.5 0.5]);
% %
% subplot(1,2,2);
% hold on; plot(Xce,FluxLim(:,end)); grid on; box on;
% set(gca,'xtick',-0.5:0.1:0.5);
% xlim([-0.5 0.5]);

%
%
%

% f4=figure(4); set(f4,'position',[1030 0 1100 420]);
% %
% subplot(1,2,1);
% hold on; plot(Xce,FluxL(:,end)); grid on;
% set(gca,'xtick',-0.5:0.1:0.5); box on;
% title('left going Flux');
% xlim([-0.5 0.5]);
% %
% subplot(1,2,2);
% hold on; plot(Xce,FluxR(:,end)); grid on; box on;
% set(gca,'xtick',-0.5:0.1:0.5);
% title('right going Flux');
% xlim([-0.5 0.5]);

end

