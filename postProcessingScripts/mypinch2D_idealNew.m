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
colormap('jet');

numProcs = 2;
filePath = '../physicsMods/pinch2D/';

numProcs = 10;
filePath = '../../fromQuartz/pinch2D_test1/';

plotBackIndex = 8; % plot time will be end-plotBackIndex

tout = loadData(filePath,numProcs,'tout');
Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Zcc = loadData(filePath,numProcs,'Zcc');
Zce = loadData(filePath,numProcs,'Zce');
rcc = loadData(filePath,numProcs,'rcc');

N  = loadData(filePath,numProcs,'N');
deltaP  = loadData(filePath,numProcs,'deltaP');
%P0  = hdf5read(thisFile,'P0');
Mx  = loadData(filePath,numProcs,'Mx');
Mz  = loadData(filePath,numProcs,'Mz');
S  = loadData(filePath,numProcs,'S');
By = loadData(filePath,numProcs,'By');
P  = loadData(filePath,numProcs,'P');
P0  = loadData(filePath,numProcs,'P0');

Vx  = loadData(filePath,numProcs,'Vx');
Vz  = loadData(filePath,numProcs,'Vz');
Jz  = loadData(filePath,numProcs,'Jz');
Ez  = loadData(filePath,numProcs,'Ez');
Ex  = loadData(filePath,numProcs,'Ex');
Cs  = loadData(filePath,numProcs,'Cs');
gamma0 = loadData(filePath,numProcs,'gamma0');
dX = Xcc(2)-Xcc(1);
dZ = Zcc(2)-Zcc(1);
Fx = loadData(filePath,numProcs,'Fx'); % - d(P+B^2/2)/dr - B^2/r

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


%%%   plot contours
%
close(figure(13));
f13=figure(13); set(f13,'position',[782 1 660 804]);
subplot(2,1,1);
pcolor(Zcc,Xcc,N(:,:,1)'); colorbar; box on
xlabel('x direction'); shading flat;
ylabel('z direction');
title('density at t = 0 t_A');
axis('equal'); axis([-0.5 0.5 0 1]);
%
%
%
subplot(2,1,2);
pcolor(Zcc,Xcc,N(:,:,end-plotBackIndex)'); colorbar; box on
%pcolor(Zcc,Xcc,P(:,:,end-plotBackIndex)'./P0'-1); colorbar; box on
xlabel('z/r_0'); shading flat;
ylabel('r/r_0'); 
title(['density at t = ',num2str(tout(end-plotBackIndex),2),' t_A']);
axis('equal'); axis([-0.5 0.5 0 1]); 
colormap('jet');


%%%   calculate global values for conservation
%
dV = rcc*2*pi*dX*dZ;
Mass = zeros(size(tout));
Entropy = zeros(size(tout));
zMomentum = zeros(size(tout));
ByFlux = zeros(size(tout));

for n=1:length(tout)
    NdV  = N(:,:,n).*dV;
    SdV  = S(:,:,n).*dV;
    MzdV = Mz(:,:,n).*dV;
    Mass(n) = sum(sum(NdV(3:end-2,3:end-2)));
    Entropy(n) = sum(sum(SdV(3:end-2,3:end-2)));
    ByFlux(n) = sum(sum(By(3:end-2,3:end-2,n)))*dX*dZ;
    zMomentum(n) = sum(sum(MzdV(3:end-2,3:end-2)));
end

close(figure(4));
figure(4); plot(tout,Mass,'displayName','mass');
hold on; plot(tout,Entropy,'displayName','entropy');
hold on; plot(tout,zMomentum,'displayName','z-momentum');
hold on; plot(tout, ByFlux,'displayName','B_y-flux');
xlabel('t/t_0'); title('global conservations');
lg4 = legend('show');
set(lg4,'location','best');




