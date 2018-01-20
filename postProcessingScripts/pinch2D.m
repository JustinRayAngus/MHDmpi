%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 2D MHD pinch simulations
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
filePath = '../physicsMods/pinch2D/';

plotBackIndex = 1; % plot time will be end-plotBackIndex

for i=1:numProcs
fileName = ['output',num2str(i-1),'.h5'];
thisFile = [filePath,fileName];
procID  = hdf5read(thisFile,'procID');
fileinfo = hdf5info(thisFile);
%
Xcc = hdf5read(thisFile,'Xcc');
Xce = hdf5read(thisFile,'Xce');
%
Zcc = hdf5read(thisFile,'Zcc');
Zce = hdf5read(thisFile,'Zce');

matW = hdf5read(thisFile,'matW');
matrixA = hdf5read(thisFile,'matrixA');
matrixB = hdf5read(thisFile,'matrixB');
matrixC = hdf5read(thisFile,'matrixC');

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
eta  = hdf5read(thisFile,'eta');
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



figure(1); hold on; pcolor(Zcc,Xcc,matW'); shading flat;
xlabel('Z'); ylabel('X'); colorbar; grid on;

% figure(4); hold on; h=pcolor(Zce,Xcc,matrixC'); shading flat;
% xlabel('Z'); ylabel('X'); colorbar; grid on;
% 
figure(2); hold on; plot(Xcc,matW(round(end/2),:));
% hold on; plot(Xcc,matrixC(round(end/2),:),'Marker','x');
%figure(2); plot(Zcc,matrixA(:,round(end/2)));
%hold on; plot(Zce,matrixC(:,round(end/2)),'Marker','x');
%display(size(matrixC));
%hold on; plot(Xcc,Xcc,'r*');

T = P./N;
Ptot = 3/2*P+B.^2/2+M.*V/2.0;
% %figure(2); hold on; plot(Xce,FluxR(:,5),'r--');
% %figure(4); hold on; plot(Xce,FluxLim(:,end),'b');
% %figure(4); hold on; plot(Xce,FluxRatio(:,end),'r');
% %%%
% figure(10); hold on; plot(Xcc,Ptot(:,end),'black','displayName','total');
% hold on; plot(Xcc,B(:,end).^2/2.0,'b','displayName','magnetic');
% hold on; plot(Xcc,1.5*P(:,end),'r','displayName','thermal');
% hold on; plot(Xcc,M(:,end).*V(:,end)/2.0,'color',[.47 .67 .19],'displayName','mean');
% title('total energy');
% box on;
if(i==1) 
    lg10=legend('show'); set(lg10,'location','best');
end
plots=0;
if(plots)
f1=figure(11); 
set(f1,'position',[1030 425 1300 840]);
%set(f1,'position',[341 436 900 840]);

subplot(2,3,1);
hold on; plot(Xcc,N(:,1),'black'); box on;
hold on; plot(Xcc,N(:,round(end/2)),'b');
hold on; plot(Xcc,N(:,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('N');
title('mass density'); axis('square');
xlim([0 1]);
%
subplot(2,3,2);
hold on; plot(Xcc,V(:,1),'black'); box on;
hold on; plot(Xcc,V(:,round(end/2)),'b');
hold on; plot(Xcc,V(:,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
xlim([0 1]);
%
subplot(2,3,3);
hold on; plot(Xcc,P(:,1),'black'); box on;
hold on; plot(Xcc,P(:,round(end/2)),'b');
hold on; plot(Xcc,P(:,end-plotBackIndex),'r'); grid on;
hold on; plot(Xcc,T(:,end-plotBackIndex),'g');
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('P');
title('thermal pressure'); axis('square');
xlim([0 1]);
%
%
subplot(2,3,4);
hold on; plot(Xce,Ez(:,1),'black'); box on;
hold on; plot(Xce,Ez(:,round(end/2)),'b');
hold on; plot(Xce,Ez(:,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Ez');
title('electric field'); axis('square');
xlim([0 1]);
%
subplot(2,3,5);
hold on; plot(Xcc,B(:,1),'black'); box on;
hold on; plot(Xcc,B(:,round(end/2)),'b'); grid on;
hold on; plot(Xcc,B(:,end-plotBackIndex),'r');
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field');
set(gca,'xtick',0.0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1]);
%
subplot(2,3,6);
hold on; plot(Xce,J(:,1),'black'); box on;
hold on; plot(Xce,J(:,round(end/2)),'b'); grid on;
hold on; plot(Xce,J(:,end-plotBackIndex),'r');
hold on; plot(Xce,J0(:,end-plotBackIndex),'g--');
xlabel('x'); ylabel('J'); axis('square');
title('current density');
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1]);
end


end

