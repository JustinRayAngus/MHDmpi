%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 2D pinch module of m=0 mode
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

plotBackIndex = 0; % plot time will be end-plotBackIndex

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
Mx  = hdf5read(thisFile,'Mx');
Mz  = hdf5read(thisFile,'Mz');
S  = hdf5read(thisFile,'S');
By  = hdf5read(thisFile,'By');
P  = hdf5read(thisFile,'P');

Vx  = hdf5read(thisFile,'Vx');
Vz  = hdf5read(thisFile,'Vz');
Jz  = hdf5read(thisFile,'Jz');
rcc = hdf5read(thisFile,'rcc');
Ez  = hdf5read(thisFile,'Ez');
Ex  = hdf5read(thisFile,'Ex');
Cs  = hdf5read(thisFile,'Cs');
gamma0 = hdf5read(thisFile,'gamma0');
FluxRatio_x  = hdf5read(thisFile,'FluxRatio_x');
FluxLim_x    = hdf5read(thisFile,'FluxLim_x');
FluxL_x    = hdf5read(thisFile,'FluxL_x');
FluxR_x    = hdf5read(thisFile,'FluxR_x');
FluxN_x  = hdf5read(thisFile,'FluxN_x');
FluxMx_x  = hdf5read(thisFile,'FluxMx_x');
FluxMz_x  = hdf5read(thisFile,'FluxMz_x');
FluxS_x  = hdf5read(thisFile,'FluxS_x');
FluxBy_x  = hdf5read(thisFile,'FluxBy_x');
tout= hdf5read(thisFile,'tout');
dX = Xcc(2)-Xcc(1);

Fx = hdf5read(thisFile,'Fx'); % - d(P+B^2/2)/dr - B^2/r

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


figure(10); hold on; plot(Xcc(3:end-2),Fx(index0,3:end-2,1)); grid on;
hold on; plot(Xcc(3:end-2),Fx(index0,3:end-2,2));
%axis([0 1 -0.002 0.002]); grid on;


plots=1;
if(plots)
f1=figure(1); 
set(f1,'position',[1030 425 1300 840]);
%set(f1,'position',[341 436 900 840]);

subplot(2,3,1);
hold on; plot(Xcc,N(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,N(index0,index1,11),'b');
hold on; plot(Xcc,N(index0,index1,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('N');
title('mass density'); axis('square');
xlim([0 1]);
%
subplot(2,3,2);
hold on; plot(Xcc,Vx(index0,index1,1),'black'); box on;
hold on; plot(Xcc,Vx(index0,index1,end/2),'b');
hold on; plot(Xcc,Vx(index0,index1,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
xlim([0 1]);
%
subplot(2,3,3);
hold on; plot(Xcc,P(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,T(index0,index1,11),'b');
hold on; plot(Xcc,P(index0,index1,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('P');
title('pressure'); axis('square');
xlim([0 1]);
ylim([0 1.2]);
%
%
subplot(2,3,4);
hold on; plot(Xcc,Ez(index0,:,1),'black'); box on;
%hold on; plot(Xce,Ez(index0,:,11),'b');
hold on; plot(Xcc,Ez(index0,:,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Ez');
title('electric field'); axis('square');
xlim([0 1]);
%
subplot(2,3,5);
hold on; plot(Xcc,By(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,By(index0,index1,11),'b'); grid on;
hold on; plot(Xcc,By(index0,index1,end-plotBackIndex),'r');
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field'); grid on;
set(gca,'xtick',0.0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1+2*dX]);
%
subplot(2,3,6);
hold on; plot(Xcc,Jz(index0,:,1),'black'); box on; grid on;
%hold on; plot(Xce,Jz(index0,:,11),'b'); grid on;
hold on; plot(Xcc,Jz(index0,:,end-plotBackIndex),'r');
xlabel('x'); ylabel('J_z'); axis('square');
title('current density');
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1]);
end

% %%%   caluculate all terms in x-momentum equation
% %
% fMx3 = -squeeze(Jzcc(3,:,:).*By(3,:,:));
% fMx2 = zeros(size(fMx3));
% fMx1 = zeros(size(fMx3));
% for k=2:length(Xce)
%     fMx1(k,:) = -(FluxMx_x(3,k,:)-FluxMx_x(3,k-1,:))/rcc(3,k)/dX;
%    % fMx2(k,:) = P(3,k,:)/rcc(3,k);
%     fMx2(k,:) = (P(3,k+1,:)+P(3,k-1,:))/2.0/rcc(3,k);
% end
% fMxtot = fMx1+fMx2+fMx3;
% % 
% figure(10); hold on; 
% plot(Xcc(3:end-2),fMxtot(3:end-2,end),'r--');
% %hold on; plot(Xcc,fMx2);
% %hold on; plot(Xcc,fMx1);


end

