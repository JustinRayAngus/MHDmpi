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
%colormap('jet');
colormap('hot');

numProcs = 2;
filePath = '../physicsMods/pinch2D/';

numProcs = 10;
%filePath = '../../fromQuartz/pinch2D/data_200nx_kR10/';
%filePath = '../../fromQuartz/pinch2D/data_400nx_kR10/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/data_800nx_kR10/'; numProcs = 40;
%filePath = '../../fromQuartz/pinch2D/data_1000nx_kR10/'; numProcs = 100;
%filePath = '../../fromQuartz/pinch2D/data_200nx_kR10_ideal1/';


%%%   kR = 1
%
filePath = '../../fromQuartz/pinch2D/kR1/data_80nz_200nx/'; numProcs = 10;
filePath = '../../fromQuartz/pinch2D/kR1/data_80nz_400nx/'; numProcs = 20;
filePath = '../../fromQuartz/pinch2D/kR1/data_160nz_800nx/'; numProcs = 80;
filePath = '../../fromQuartz/pinch2D/kR1/data_320nz_800nx/'; numProcs = 80;

%%%   kR = 5
%
%filePath = '../../fromQuartz/pinch2D/kR5/data_80nz_200nx/'; numProcs = 10;
%filePath = '../../fromQuartz/pinch2D/kR5/data_80nz_400nx/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/kR5/data_160nz_800nx/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR5/data_320nz_800nx/'; numProcs = 80;


%%%   kR = 10
%
%filePath = '../../fromQuartz/pinch2D/kR10/data_80nz_400nx/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/kR10/data_160nz_800nx/'; numProcs = 80;
filePath = '../../fromQuartz/pinch2D/kR10/data_400nz_800nx/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/testing/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/kR10/testing2/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/testing3/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/testing4/'; numProcs = 80;


plotBackIndex = 50; % plot time will be end-plotBackIndex
thist = 9.4; %7.9;
thist2 = 10.2; %9.4;

t0 = 3.6137e-8; % see normalizationParameters_zPinch.m
tA = 2.5e-8;    % from PoP 17, 072107 (2010)
%tA = 5.4036e-08;  % r0/VTi
%tA = 3.8209e-08;


tout = loadData(filePath,numProcs,'tout');
tout = tout*t0/tA; % normalize to 2010 paper alfven time
[~,tindex] = min(abs(tout-thist));
[~,tindex2] = min(abs(tout-thist2));

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


% figure(10); hold on; plot(Xcc(3:end-2),Fx(index0,3:end-2,1)); grid on;
% hold on; plot(Xcc(3:end-2),Fx(index0,3:end-2,2));
% %axis([0 1 -0.002 0.002]); grid on;


plots=1;
if(plots)
f1=figure(1); 
set(f1,'position',[1030 425 1300 840]);
%set(f1,'position',[341 436 900 840]);

subplot(2,3,1);
hold on; plot(Xcc,N(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,N(index0,index1,11),'b');
hold on; plot(Xcc,N(index0,index1,tindex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('N');
title('mass density'); axis('square');
xlim([0 1]);
%
subplot(2,3,2);
hold on; plot(Xcc,Vx(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,Vx(index0,index1,end/2),'b');
hold on; plot(Xcc,Vx(index0,index1,tindex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
xlim([0 1]);
%
subplot(2,3,3);
hold on; plot(Xcc,P(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,T(index0,index1,11),'b');
hold on; plot(Xcc,P(index0,index1,tindex),'r'); grid on;
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
hold on; plot(Xcc,Ez(index0,:,tindex),'r'); grid on;
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Ez');
title('electric field'); axis('square');
xlim([0 1]);
%
subplot(2,3,5);
hold on; plot(Xcc,By(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,By(index0,index1,11),'b'); grid on;
hold on; plot(Xcc,By(index0,index1,tindex),'r');
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field'); grid on;
set(gca,'xtick',0.0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1+2*dX]);
%
subplot(2,3,6);
hold on; plot(Xcc,Jz(index0,:,1),'black'); box on; grid on;
%hold on; plot(Xce,Jz(index0,:,11),'b'); grid on;
hold on; plot(Xcc,Jz(index0,:,tindex),'r');
xlabel('x'); ylabel('J_z'); axis('square');
title('current density');
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1]);
end


%%%   plot contours
%
close(figure(15));
f13=figure(15); set(f13,'position',[1282 300 660 804]);
subplot(2,1,1);
pcolor(Zcc,Xcc,N(:,:,tindex)'); colorbar; box on
xlabel('x direction'); shading flat;
ylabel('z direction');
title(['density at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
subplot(2,1,2);
pcolor(Zcc,Xcc,N(:,:,tindex2)'); colorbar; box on
%pcolor(Zcc,Xcc,P(:,:,end-plotBackIndex)'./P0'-1); colorbar; box on
xlabel('z/r_0'); shading flat;
ylabel('r/r_0'); 
title(['density at t = ',num2str(tout(tindex2),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
colormap('jet');


%%%   calculate global values for conservation
%
dV = rcc*2*pi*dX*dZ;
Mass = zeros(size(tout));
Entropy = zeros(size(tout));
zMomentum = zeros(size(tout));
ByFlux = zeros(size(tout));
Energy = zeros(size(tout));
Edens = 0.5*(Mx.*Vx+Mz.*Vz) + 1.5*P + By.^2/2.0; 

for n=1:length(tout)
    NdV  = N(:,:,n).*dV;
    SdV  = S(:,:,n).*dV;
    MzdV = Mz(:,:,n).*dV;
    EdensdV = Edens.*dV;
    Mass(n) = sum(sum(NdV(3:end-2,3:end-2)));
    Entropy(n) = sum(sum(SdV(3:end-2,3:end-2)));
    ByFlux(n) = sum(sum(By(3:end-2,3:end-2,n)))*dX*dZ;
    zMomentum(n) = sum(sum(MzdV(3:end-2,3:end-2)));
    Energy(n) = sum(sum(EdensdV(3:end-2,3:end-2)));
end


close(figure(4));
figure(4); plot(tout,Mass,'displayName','mass');
hold on; plot(tout,Entropy,'displayName','entropy');
hold on; plot(tout,zMomentum,'displayName','z-momentum');
hold on; plot(tout, ByFlux,'displayName','B_y-flux');
xlabel('t/t_0'); title('global conservations');
lg4 = legend('show');
set(lg4,'location','best');


%%%   calculate perturbed mass density per unit length 
%
dA = rcc*2*pi*dX;
deltaN = zeros(size(N));
Dz = zeros(length(Zcc)-4,length(tout)); % mass per unit length z


for n=1:length(tout)
  %  deltaN(:,:,n) = N(:,:,n);
    deltaN(:,:,n) = S(:,:,n);
    deltaNdA = squeeze(deltaN(3:end-2,3:end-2,n).*dA(3:end-2,3:end-2));
    Dz(:,n) = sum(deltaNdA,2);
end
%Dz = squeeze(N(3:end-2,60,:));


%%%   calculate Fourier transform of all modes and amplitudes
%%%   assumed period on Zce grid
%
Nz = length(Zce)-2;
L = Zce(end-1)-Zce(2);

Zmodes = 0:1:(Nz-1)/2; % assuming N is odd
k = 2*pi*Zmodes/L;
DzFT = zeros(length(k),length(tout)); % Fourier transform

for n=1:length(tout)
   for j=1:length(Zmodes)
      Dzexp = Dz(:,n).*exp(-1i*k(j)*Zcc(3:end-2));
      DzFT(j,n) = sum(Dzexp)*dZ;
   end
end
DzFT_amp = abs(DzFT)/Mass(1);
%DzFT_amp = abs(real(DzFT))/Mass(1);
DzFT_amp_odd = abs(imag(DzFT))/Mass(1);

%tout = sqrt(2)*tout;
f11=figure(14); set(f11,'position',[1000 200 510 850]);

subplot(2,1,1);
semilogy(tout,DzFT_amp(2,:),'displayName',['kR=',num2str(k(2))]);
hold on; plot(tout,DzFT_amp(3,:),'displayName',['kR=',num2str(k(3))]);
hold on; plot(tout,DzFT_amp(4,:),'displayName',['kR=',num2str(k(4))]);
xlabel('t/t_A'); ylabel('mode amplitude'); 
title('Fourier Amplitudes'); grid on, grid off; grid on;
lg9=legend('show'); set(lg9,'location','southeast');
ylim([1.0e-10 1]); axis('square');
set(gca,'ytick',10.^(-10:2:0));

subplot(2,1,2);
semilogy(tout,DzFT_amp_odd(2,:),'displayName',['kR=',num2str(k(2))]);
hold on; plot(tout,DzFT_amp_odd(3,:),'displayName',['kR=',num2str(k(3))]);
hold on; plot(tout,DzFT_amp_odd(4,:),'displayName',['kR=',num2str(k(4))]);
xlabel('t/t_A'); ylabel('mode amplitude'); 
title('odd Fourier Amplitudes'); grid on, grid off; grid on;
lg9=legend('show'); set(lg9,'location','southeast');
ylim([1.0e-10 1]); axis('square');
set(gca,'ytick',10.^(-10:2:0));

%%% for looking at modes of N
%hold on; plot(tout,2.8e-5*10.^(tout/2.15),'b--');
%hold on; plot(tout,2.0e-10*10.^(tout*2/2.15),'r--');
%hold on; plot(tout,7.5e-19*10.^(tout*3/2.15),'linestyle','--'); ? should
%be 3 times as fast, not 4 times as fast


%%% for looking at modes of S
%
%hold on; plot(tout,4.2e-5*10.^(tout/2.15),'b--');
%hold on; plot(tout,8.5e-10*10.^(tout*2/2.15),'r--');
%hold on; plot(tout,2.0e-14*10.^(tout*3/2.15),'linestyle','--');


%%%   calculate slope of Fourier modes
%
mDzFT_amp = zeros(length(k),length(tout)-1);
for i=1:length(k)
    for j=1:length(tout)-1
        mDzFT_amp(i,j) = (log(DzFT_amp(i,j+1))-log(DzFT_amp(i,j)))/(tout(j+1)-tout(j));
    end
end


figure(23); hold on; plot(tout(1:end-1),mDzFT_amp(2,:)); box on;
hold on; plot(tout(1:end-1),mDzFT_amp(3,:));
hold on; plot(tout(1:end-1),mDzFT_amp(4,:));
xlabel('t/t_A'); ylabel('growth rate');
title('derivative of log Fourier amps');
axis([0 10 -1 5]); grid on;



NperZ = Mass(1)/dZ/length(Zcc(3:end-2))/pi; % ave density per unit length
Mi = 1.008*1.6605402e-27; % [kg]
n0 = 1.0e24; % [1/m^3]
R0 = 5.0e-3; % [m]

rhobar = Mi*NperZ*n0; % [kg/m^3]





