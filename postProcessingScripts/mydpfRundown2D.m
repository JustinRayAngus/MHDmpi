%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 2D dpfRundown module
%%%
%%%   see G.A. Sod,. J. Comp. Phys. 27, 1-31 (1978)
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
filePath = '../physicsMods/dpfRundown2D/';

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
Jz0  = hdf5read(thisFile,'Jz0'); % curl(B)
Ez  = hdf5read(thisFile,'Ez');
Ex  = hdf5read(thisFile,'Ex');
Cs  = hdf5read(thisFile,'Cs');
eta  = hdf5read(thisFile,'eta');
gamma0 = hdf5read(thisFile,'gamma0');
FluxRatio_x  = hdf5read(thisFile,'FluxRatio_x');
FluxLim_x    = hdf5read(thisFile,'FluxLim_x');
FluxL_x    = hdf5read(thisFile,'FluxL_x');
FluxR_x    = hdf5read(thisFile,'FluxR_x');
FluxN_x  = hdf5read(thisFile,'FluxN_x');
FluxMx_x  = hdf5read(thisFile,'FluxMx_x');
FluxS_x  = hdf5read(thisFile,'FluxS_x');
FluxBy_x  = hdf5read(thisFile,'FluxBy_x');
FluxEz_x  = hdf5read(thisFile,'FluxEz_x');
tout= hdf5read(thisFile,'tout');
dX = Xcc(2)-Xcc(1);


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


figure(10); hold on; plot(Xcc,Ptot(index0,index1,end),'black','displayName','total');
hold on; plot(Xcc,By(index0,index1,end).^2/2.0,'b','displayName','magnetic');
hold on; plot(Xcc,P(index0,index1,end),'r','displayName','thermal');
hold on; plot(Xcc,Mx(index0,index1,end).*Vx(index0,index1,end)/2.0,'color',[.47 .67 .19],'displayName','mean');
title('total energy');
box on;
if(i==1) 
    lg10=legend('show'); set(lg10,'location','best');
end
plots=1;
if(plots)
f1=figure(11); 
set(f1,'position',[1030 425 1300 840]);
%set(f1,'position',[341 436 900 840]);

subplot(2,3,1);
hold on; plot(Xcc,N(index0,index1,1),'black'); box on;
hold on; plot(Xcc,N(index0,index1,11),'b');
hold on; plot(Xcc,N(index0,index1,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('N');
title('mass density'); axis('square');
xlim([0 1]);
%
subplot(2,3,2);
hold on; plot(Xcc,Vx(index0,index1,1),'black'); box on;
hold on; plot(Xcc,Vx(index0,index1,11),'b');
hold on; plot(Xcc,Vx(index0,index1,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
xlim([0 1]);
%
subplot(2,3,3);
hold on; plot(Xcc,P(index0,index1,1),'black'); box on;
hold on; plot(Xcc,P(index0,index1,11),'b');
hold on; plot(Xcc,P(index0,index1,end-plotBackIndex),'r'); grid on;
hold on; plot(Xcc,T(index0,index1,end-plotBackIndex),'g');
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('P');
title('thermal pressure'); axis('square');
xlim([0 1]);
%
%
subplot(2,3,4);
hold on; plot(Xce,Ez(index0,:,1),'black'); box on;
hold on; plot(Xce,Ez(index0,:,11),'b');
hold on; plot(Xce,Ez(index0,:,end-plotBackIndex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Ez');
title('electric field'); axis('square');
xlim([0 1]);
%
subplot(2,3,5);
hold on; plot(Xcc,By(index0,index1,1),'black'); box on;
hold on; plot(Xcc,By(index0,index1,11),'b'); grid on;
hold on; plot(Xcc,By(index0,index1,end-plotBackIndex),'r');
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field');
set(gca,'xtick',0.0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1]);
%
subplot(2,3,6);
hold on; plot(Xce,Jz(index0,:,1),'black'); box on;
hold on; plot(Xce,Jz(index0,:,11),'b'); grid on;
hold on; plot(Xce,Jz(index0,:,end-plotBackIndex),'r');
hold on; plot(Xce,Jz0(index0,:,end-plotBackIndex),'g--');
xlabel('x'); ylabel('J'); axis('square');
title('current density');
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1]);
end


end

