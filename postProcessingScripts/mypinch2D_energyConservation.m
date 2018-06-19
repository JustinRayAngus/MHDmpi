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
newDeck=0;



filePath = '../../fromQuartz/pinch2D/entropy_v0/Li0.0/ka3.0/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v0/Li1.5e-2/kR3.0_stable/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v0/Li1.5e-2/ka30_Te/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v0/Li0.12/ka3.0/'; numProcs = 20; newDeck = 2;


filePath = '../../fromQuartz/pinch2D/driftIdeal_v2/ka3.0_1/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/driftIdeal_v2/ka3.0_2/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/driftIdeal_v2/ka3.0_3/'; numProcs = 20; newDeck = 2;


plotBackIndex = 1; % plot time will be end-plotBackIndex
thist = 15; %9.4; %7.9;
thist2 = 20; %10.2; %9.4;

%t0 = 3.6137e-8; % see normalizationParameters_zPinch.m
t0 = 1.2046e-8; % using pinch radius for length scale
tA = 2.5e-8;    % from PoP 17, 072107 (2010)
%tA = 5.4036e-08;  % r0/VTi
%tA = 3.8209e-08;


tout = loadData(filePath,numProcs,'tout');
%tout = tout*t0/tA; % normalize to 2010 paper alfven time
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
      Ee_source = loadData(filePath,numProcs,'Ee_source');
      Ei_source = loadData(filePath,numProcs,'Ei_source');
      Te = loadData(filePath,numProcs,'Te');
      Ti = loadData(filePath,numProcs,'Ti');
      Vex = Vx - lambda0*Jx./N;
      Vez = Vz - lambda0*Jz./N;
      Ee  = loadData(filePath,numProcs,'Ee');
      Ei  = loadData(filePath,numProcs,'Ei');
  end
end


%%%   plot contours
%
close(figure(12));
f13=figure(12); set(f13,'position',[1282 300 660 804]);
subplot(2,1,1);
pcolor(Zcc,Xcc,N(:,:,tindex)'); colorbar; box on
xlabel('x direction'); shading flat;
ylabel('z direction');
title(['density at t = ',num2str(tout(tindex),3),' t_A']);
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
Energy_By = zeros(size(tout));
Energy_Pe = zeros(size(tout));
Energy_Pi = zeros(size(tout));
Energy_Ui = zeros(size(tout));
Edens = 0.5*(Mx.*Vx + Mz.*Vz) + P/(gamma0-1) + By.^2/2.0; 
%
Ee0 = zeros(size(tout));
Ei0 = zeros(size(tout));
Se0 = zeros(size(tout));
Si0 = zeros(size(tout));



for n=1:length(tout)
    NdV  = N(:,:,n).*dV;
    SdV  = P(:,:,n)./N(:,:,n).^gamma0.*dV;
    MzdV = Mz(:,:,n).*dV;
    EdensdV_By = By(:,:,n).^2/2.0.*dV;
    EdensdV_Pe = Te(:,:,n).*N(:,:,n).*dV/(gamma0-1);
    EdensdV_Pi = Ti(:,:,n).*N(:,:,n).*dV/(gamma0-1);
    EdensdV_Ui = (Mx(:,:,n).*Vx(:,:,n)+Mz(:,:,n).*Vz(:,:,n)).*dV/2;
    EdensdV = Edens(:,:,n).*dV;
    %
    EedV = Ee(:,:,n).*dV;
    EidV = Ei(:,:,n).*dV;
    SedV = Ee_source(:,:,n).*dV;
    SidV = Ei_source(:,:,n).*dV;
    %
    Mass(n) = sum(sum(NdV(3:end-2,3:end-2)));
    Entropy(n) = sum(sum(SdV(3:end-2,3:end-2)));
    ByFlux(n) = sum(sum(By(3:end-2,3:end-2,n)))*dX*dZ;
    zMomentum(n) = sum(sum(MzdV(3:end-2,3:end-2)));
    Energy(n) = sum(sum(EdensdV(3:end-2,3:end-2)));
    Energy_By(n) = sum(sum(EdensdV_By(3:end-2,3:end-2)));
    Energy_Pe(n) = sum(sum(EdensdV_Pe(3:end-2,3:end-2)));
    Energy_Pi(n) = sum(sum(EdensdV_Pi(3:end-2,3:end-2)));    
    Energy_Ui(n) = sum(sum(EdensdV_Ui(3:end-2,3:end-2)));
    %
    Ee0(n) = sum(sum(EedV(3:end-2,3:end-2)));
    Ei0(n) = sum(sum(EidV(3:end-2,3:end-2)));
    Se0(n) = sum(sum(SedV(3:end-2,3:end-2)));
    Si0(n) = sum(sum(SidV(3:end-2,3:end-2)));
end


close(figure(4));
figure(4); %plot(tout,Mass,'displayName','mass');
plot(tout, Ee0,'displayName','Ee','color','b');
hold on; plot(tout, Ei0,'displayName','Ei','color','r');
hold on; plot(tout,Entropy,'displayName','entropy');
%hold on; plot(tout,zMomentum,'displayName','z-momentum');
%hold on; plot(tout, ByFlux,'displayName','B_y-flux');
hold on; plot(tout, Energy,'displayName','Energy');
xlabel('t/t_0'); title('global conservations');
lg4 = legend('show');
set(lg4,'location','southwest');



figure(11); 
hold on; plot(tout,Energy_Ui); box on;
hold on; plot(tout,Energy_Pi-Energy_Pi(1));
hold on; plot(tout,Energy_Pe-Energy_Pe(1));
hold on; plot(tout,Energy_By-Energy_By(1));
title('\Delta Energy partition');
lg11=legend('ion mean','ion therm','ele therm','magnetic');
set(lg11,'location','best');

figure(5); plot(tout,Se0,'b');
hold on; plot(tout,Si0,'r');



%%%   calculate expected change in energy and compare with actual change
%
Se0tot = cumtrapz(tout,Se0);
Si0tot = cumtrapz(tout,Si0);


close(figure(7));
f7=figure(7); set(f7,'position',[1400 50 1000 420]);
subplot(1,2,1);
plot(tout,Ei0-Ei0(1));
hold on; plot(tout,Si0tot,'r--');
lg7=legend('Ei0-Ei0(1)','\int Si0 dt');
set(lg7,'location','best');
title('ion energy: MiU^2/2 + 3Pi/2');

subplot(1,2,2);
plot(tout,Ee0-Ee0(1));
hold on; plot(tout,Se0tot,'r--');
lg8=legend('Ee0-Ee0(1)','\int Se0 dt');
set(lg8,'location','best');
title('ele+B energy: 3Pe/2 + B^2/2');



