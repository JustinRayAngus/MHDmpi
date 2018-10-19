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
%colormap('hot');

numProcs = 2;
filePath = '../physicsMods/pinch2D/';
newDeck=0;

numProcs = 10;
%filePath = '../../fromQuartz/pinch2D/data_200nx_kR10/';
%filePath = '../../fromQuartz/pinch2D/data_400nx_kR10/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/data_800nx_kR10/'; numProcs = 40;
%filePath = '../../fromQuartz/pinch2D/data_1000nx_kR10/'; numProcs = 100;
%filePath = '../../fromQuartz/pinch2D/data_200nx_kR10_ideal1/';


%%%   kR = 1
%
%filePath = '../../fromQuartz/pinch2D/kR1/data_80nz_200nx/'; numProcs = 10;
%filePath = '../../fromQuartz/pinch2D/kR1/data_80nz_400nx/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/kR1/data_160nz_800nx/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR1/data_320nz_800nx/'; numProcs = 80;

%%%   kR = 5
%
filePath = '../../fromQuartz/pinch2D/kR5/data_80nz_200nx/'; numProcs = 10;
filePath = '../../fromQuartz/pinch2D/kR5/data_80nz_400nx/'; numProcs = 20;
filePath = '../../fromQuartz/pinch2D/kR5/data_160nz_800nx/'; numProcs = 80;
filePath = '../../fromQuartz/pinch2D/kR5/data_320nz_800nx/'; numProcs = 80;


%%%   kR = 20
%
%filePath = '../../fromQuartz/pinch2D/kR20/data_80nz_400nx/'; numProcs = 20;


%%%   kR = 30
%
%filePath = '../../fromQuartz/pinch2D/kR30/data_80nz_400nx/'; numProcs = 20;


%%%   kR = 40
%
%filePath = '../../fromQuartz/pinch2D/kR40/data_80nz_400nx/'; numProcs = 20;


%%%   kR = 10
%
%filePath = '../../fromQuartz/pinch2D/kR10/data_80nz_400nx/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/kR10/data_160nz_800nx/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/data_400nz_800nx/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/testing/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/kR10/testing2/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/testing3/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/testing4/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/data_80nz_400nx_ideal1/'; numProcs = 20;
%
% ideal2/testing1-3 is with epsilon = 1.0e-4
%
%filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing1/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing2/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing3/'; numProcs = 20; newDeck = 1;
%
% ideal2/testing4-6 is with epsilon = 1.0e-6
%
%filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing4/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing5/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing6/'; numProcs = 20; newDeck = 1;
%
% Hall_v1/
%
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test4_Li=0/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test4_Li=1.0e-3/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test5_Li=0/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test5_Li=1.0e-3/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test5_Li=4.6e-3/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test5_Li=1.0e-2/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test6_Li=0/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test6_Li=1.0e-3/'; numProcs = 20; newDeck = 1;
%
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v1/test4_Li=1.0e-3/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/testing6/'; numProcs = 20; newDeck = 1;

%
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_tauei/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable_tauei/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-2/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=2.0e-2_again/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v0/kR40/test5_Li=5.0e-3/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v0/kR40/test5_Li=5.0e-3_nuT/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v0/Li0.0/ka100/'; numProcs = 20; newDeck = 2;


%filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.17/ka3.0_taui1.0e-3/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/ka3.0_noGyroVisc_taui1.0e-2_smallerDt/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/ka3.0_noGyroVisc_taui1.0e-2/oldData/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/ka3.0_noGyroVisc/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/noGyroVisc/ka3.0_taui1.0e-2/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/noGyroVisc/ka3.0_taui1.0e-2_nuTherm100/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/noGyroVisc/ka0.3_taui1.0e-2/'; numProcs = 20; newDeck = 2;
%
filePath = '/Volumes/LaCie/zPinch/myMHD/entropy_v1/Li0.12/noGyroVisc/ka3.0_taui1.0e-2/'; numProcs = 20; newDeck = 2;
%filePath = '/Volumes/LaCie/zPinch/myMHD/entropy_v1/Li0.12/noGyroVisc/ka3.0_taui1.0e-2_nuTherm1/'; numProcs = 20; newDeck = 2;
%filePath = '/Volumes/LaCie/zPinch/myMHD/entropy_v1/Li0.12/noGyroVisc/ka3.0_testing/'; numProcs = 20; newDeck = 2;

%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/ka3.0_HallC2/'; numProcs = 20; newDeck = 2;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/ka3.0_HallTVD/'; numProcs = 20; newDeck = 2;
filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/ka3.0_taui1.0e-2/'; numProcs = 20; newDeck = 2;

filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/ka10.0_nuTherm10/'; numProcs = 20; newDeck = 2;
filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/ka3.0_nuTherm0/'; numProcs = 20; newDeck = 2;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/withShear/testing3/'; numProcs = 20; newDeck = 2;
filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/withShear/ka3.0_M0.5_ideal/'; numProcs = 20; newDeck = 2;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/withShear/ka3.0_M0.5_taui1.0e-3_splitCspeed/'; numProcs = 20; newDeck = 2;

filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/Paraschiv_Fig10/kR5.0_M1.25/'; numProcs = 20; newDeck = 2;



%filePath = '../../fromQuartz/pinch2D/entropy_v1/Li1.5e-2/ka3.0_new/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v2/testing4_taui1.0e-3/'; numProcs = 20; newDeck = 2;


%filePath = '../../fromQuartz/pinch2D/driftIdeal_v2/ka3.0_0/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/driftIdeal_v2/ka3.0_0_Li0p015/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/driftIdeal_v2/ka3.0_3_Li0p015/'; numProcs = 20; newDeck = 2;



plotBackIndex = 1; % plot time will be end-plotBackIndex
thist = 10; %9.4; %7.9;
thist2 = 14.2; %10.2; %9.4;

%t0 = 3.6137e-8; % see normalizationParameters_zPinch.m
t0 = 1.2040e-8; % using pinch radius for length scale
tA = 2.5e-8;    % from PoP 17, 072107 (2010)
%tA = 5.4036e-08;  % r0/VTi
%tA = 3.8209e-08;


tout = loadData(filePath,numProcs,'tout');
tout = tout*t0/tA; % normalize to 2010 paper alfven time
[~,tindex] = min(abs(tout-thist2));
[~,tindex2] = min(abs(tout-tout(end-1)));

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
Jx  = loadData(filePath,numProcs,'Jx');
Ez  = loadData(filePath,numProcs,'Ez');
Ex  = loadData(filePath,numProcs,'Ex');
Cs  = loadData(filePath,numProcs,'Cs');
gamma0 = loadData(filePath,numProcs,'gamma0');
dX = Xcc(2)-Xcc(1);
dZ = Zcc(2)-Zcc(1);
Fx = loadData(filePath,numProcs,'Fx'); % - d(P+B^2/2)/dr - B^2/r

Jz0 = Jz;
Ez0 = Ez;
if(newDeck>=1)
  Jz0  = loadData(filePath,numProcs,'Jz0');
  Ez0  = loadData(filePath,numProcs,'Ez0');
  if(newDeck>=2)
      delta0 = loadData(filePath,numProcs,'delta0');
      lambda0 = loadData(filePath,numProcs,'lambda0');
   %   Ee_source = loadData(filePath,numProcs,'Ee_source');
   %   Ei_source = loadData(filePath,numProcs,'Ei_source');
      Te = loadData(filePath,numProcs,'Te');
      Ti = loadData(filePath,numProcs,'Ti');
      Vex = Vx - lambda0*Jx./N;
      Vez = Vz - lambda0*Jz./N;
   %   Ee  = loadData(filePath,numProcs,'Ee');
   %   Ei  = loadData(filePath,numProcs,'Ei');
  end
end


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
%set(gca,'xtick',0:0.25:3);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('N');
title('mass density'); axis('square');
xlim([0 Xce(end-1)]);
%
subplot(2,3,2);
hold on; plot(Xcc,Vx(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,Vx(index0,index1,end/2),'b');
hold on; plot(Xcc,Vx(index0,index1,tindex),'r'); grid on;
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('V');
title('velocity'); axis('square');
xlim([0 Xce(end-1)]);
%
subplot(2,3,3);
hold on; plot(Xcc,P(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,T(index0,index1,11),'b');
hold on; plot(Xcc,P(index0,index1,tindex),'r'); grid on;
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('P');
title('pressure'); axis('square');
xlim([0 Xce(end-1)]);
ylim([0 1.2]);
%
%
subplot(2,3,4);
hold on; plot(Xcc,Ez(index0,:,1),'black'); box on;
%hold on; plot(Xce,Ez(index0,:,11),'b');
hold on; plot(Xcc,Ez(index0,:,tindex),'r'); grid on;
hold on; plot(Xcc,Ez0(index0,:,tindex),'g--');
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Ez');
title('electric field'); axis('square');
xlim([0 Xce(end-1)]);
%
subplot(2,3,5);
hold on; plot(Xcc,By(index0,index1,1),'black'); box on;
%hold on; plot(Xcc,By(index0,index1,11),'b'); grid on;
hold on; plot(Xcc,By(index0,index1,tindex),'r');
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field'); grid on;
%set(gca,'xtick',0.0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 Xce(end-1)]);
%
subplot(2,3,6);
hold on; plot(Xcc,Jz(index0,:,1),'black'); box on; grid on;
%hold on; plot(Xce,Jz(index0,:,11),'b'); grid on;
hold on; plot(Xcc,Jz(index0,:,tindex),'r');
hold on; plot(Xcc,Jz0(index0,:,tindex),'g--');
xlabel('x'); ylabel('J_z'); axis('square');
title('current density');
%set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 Xce(end-1)]);
end


%%%   plot contours
%
close(figure(15));
f13=figure(15); set(f13,'position',[1282 300 660 804]);
subplot(2,1,1);
pcolor(Zcc,Xcc,N(:,:,tindex)'); colorbar; box on
xlabel('x direction'); shading flat;
ylabel('z direction');
title(['density at t = ',num2str(tout(tindex),3),' t_0']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(2) Xce(end-1)]);
%
%
%
subplot(2,1,2);
pcolor(Zcc,Xcc,N(:,:,tindex2)'); colorbar; box on
%pcolor(Zcc,Xcc,P(:,:,end-plotBackIndex)'./P0'-1); colorbar; box on
xlabel('z/r_0'); shading flat;
ylabel('r/r_0'); 
title(['density at t = ',num2str(tout(tindex2),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(2) Xce(end-1)]);
colormap('jet');


%%%   calculate global values for conservation
%
dV = rcc*2*pi*dX*dZ;
Mass = zeros(size(tout));
Entropy = zeros(size(tout));
zMomentum = zeros(size(tout));
ByFlux = zeros(size(tout));
Energy = zeros(size(tout));
%Energy2 = zeros(size(tout));
Edens = 0.5*(Mx.*Vx + Mz.*Vz) + P/(gamma0-1) + By.^2/2.0; 
%Edens2 = Ee + Ei;
S = N.*(Te+Ti)./N.^gamma0;

for n=1:length(tout)
    NdV  = N(:,:,n).*dV;
    SdV  = N(:,:,n).*log(S(:,:,n)).*dV;
    MzdV = Mz(:,:,n).*dV;
    EdensdV = Edens(:,:,n).*dV;
    %Edens2dV = Edens2(:,:,n).*dV;
    Mass(n) = sum(sum(NdV(3:end-2,3:end-2)));
    Entropy(n) = sum(sum(SdV(3:end-2,3:end-2)));
    ByFlux(n) = sum(sum(By(3:end-2,3:end-2,n)))*dX*dZ;
    zMomentum(n) = sum(sum(MzdV(3:end-2,3:end-2)));
    Energy(n) = sum(sum(EdensdV(3:end-2,3:end-2)));
    %Energy2(n) = sum(sum(Edens2dV(3:end-2,3:end-2)));
end


close(figure(4));
figure(4); plot(tout,Mass,'displayName','mass');
hold on; plot(tout,Entropy,'displayName','entropy');
hold on; plot(tout,zMomentum,'displayName','z-momentum');
hold on; plot(tout, ByFlux,'displayName','B_y-flux');
hold on; plot(tout, Energy,'displayName','Energy');
%hold on; plot(tout, Energy2,'r--');
xlabel('t/t_0'); title('global conservations');
lg4 = legend('show');
set(lg4,'location','best');


%%%   calculate perturbed mass density per unit length 
%
dA = rcc*2*pi*dX;
deltaN = zeros(size(N));
Dz = zeros(length(Zcc)-4,length(tout)); % mass per unit length z


for n=1:length(tout)
    deltaN(:,:,n) = N(:,:,n);
 %   deltaN(:,:,n) = Te(:,:,n);
    %  deltaN(:,:,n) = S(:,:,n);
    deltaNdA = squeeze(deltaN(3:end-2,3:end-2,n).*dA(3:end-2,3:end-2));
    %deltaNdA = squeeze(deltaN(3:end-2,50:350,n).*dA(3:end-2,50:350));
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
f10=figure(10); set(f10,'position',[1000 200 510 850]);

subplot(2,1,1);
semilogy(tout,DzFT_amp(2,:),'displayName',['ka=',num2str(k(2))]);
hold on; plot(tout,DzFT_amp(3,:),'displayName',['ka=',num2str(k(3))]);
hold on; plot(tout,DzFT_amp(4,:),'displayName',['ka=',num2str(k(4))]);
xlabel('t/t_A'); ylabel('mode amplitude'); 
title('Fourier Amplitudes'); grid on, grid off; grid on;
lg9=legend('show'); set(lg9,'location','southeast');
ylim([1.0e-10 1]); axis('square');
set(gca,'ytick',10.^(-10:2:0));

subplot(2,1,2);
semilogy(tout,DzFT_amp_odd(2,:),'displayName',['ka=',num2str(k(2))]);
hold on; plot(tout,DzFT_amp_odd(3,:),'displayName',['ka=',num2str(k(3))]);
hold on; plot(tout,DzFT_amp_odd(4,:),'displayName',['ka=',num2str(k(4))]);
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
xlabel('t/t_0'); ylabel('growth rate');
title('derivative of log Fourier amps');
axis([0 tout(end) -1 5]); grid on;



NperZ = Mass(1)/dZ/length(Zcc(3:end-2))/pi; % ave density per unit length
Mi = 1.008*1.6605402e-27; % [kg]
n0 = 1.0e24; % [1/m^3]
R0 = 5.0e-3; % [m]

rhobar = Mi*NperZ*n0; % [kg/m^3]





