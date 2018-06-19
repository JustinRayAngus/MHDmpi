%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   something funny going on at r=0 for Jz0 = d(rB)/dr/r
%%%   It's a BC issue....
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_tauei/'; numProcs = 20; newDeck = 2;


thist = 10; %10.2; %9.4;

t0 = 3.6137e-8; % see normalizationParameters_zPinch.m
tA = 2.5e-8;    % from PoP 17, 072107 (2010)
%tA = 5.4036e-08;  % r0/VTi
%tA = 3.8209e-08;


tout = loadData(filePath,numProcs,'tout');
tout = tout*t0/tA; % normalize to 2010 paper alfven time
[~,tindex] = min(abs(tout-thist));

Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Zcc = loadData(filePath,numProcs,'Zcc');
Zce = loadData(filePath,numProcs,'Zce');
rcc = loadData(filePath,numProcs,'rcc');


By = loadData(filePath,numProcs,'By');
Jz = loadData(filePath,numProcs,'Jz');
Ez = loadData(filePath,numProcs,'Ez');
FluxEz_x = loadData(filePath,numProcs,'FluxEz_x');

dX = Xcc(2)-Xcc(1);
dZ = Zcc(2)-Zcc(1);
Jz0 = Jz;
Ez0 = Ez;
if(newDeck>=1)
  Jz0  = loadData(filePath,numProcs,'Jz0');
  Ez0  = loadData(filePath,numProcs,'Ez0');
end


index0 = round(length(Ez(:,1,1)))/2.0;

f1=figure(1); 
set(f1,'position',[1130 800 1370 500]);
%
subplot(1,3,1);
hold on; plot(Xcc,Ez(index0,:,tindex),'black'); grid on; box on;
hold on; plot(Xcc,Ez0(index0,:,tindex),'r--');
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x'); ylabel('Ez');
title('electric field'); axis('square');
xlim([0 1]);
%
subplot(1,3,2);
hold on; plot(Xcc,By(index0,:,tindex),'black'); box on;
xlabel('x'); ylabel('B'); axis('square');
title('magnetic field'); grid on;
set(gca,'xtick',0.0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1+2*dX]);
%
subplot(1,3,3);
hold on; plot(Xcc,Jz(index0,:,tindex),'black'); box on; grid on;
hold on; plot(Xcc,Jz0(index0,:,tindex),'r--');
xlabel('x'); ylabel('J_z'); axis('square');
title('current density');
set(gca,'xtick',0:0.25:2);
%set(gca,'ytick',1:0.5:3);
xlim([0 1]);


for i=1:length(Xcc)
   By(:,i,:) = Xcc(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   calculate Jz00 = d(rB)/dr/r
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%   compute cell-center flux = r*B
%
FluxEzcc = zeros(size(By));
for i=1:length(By(1,:,1))
    FluxEzcc(:,i,:) = Xcc(i)*By(:,i,:);
end
FluxEzcc(:,1,:) = 0;
FluxEzcc(:,2,:) = 0;


%%%   compute cell-edge flux = r*B
%
FluxEzce = zeros(size(FluxEz_x));
for i=1:length(Xce)
    FluxEzce(:,i,:) = (FluxEzcc(:,i+1,:) + FluxEzcc(:,i,:))/2.0;
end


%%%   plot fluxes near r=0
%
f2=figure(2); 
hold on; plot(Xcc,FluxEzcc(index0,:,tindex),'displayName','Flux_c_c');
% hold on; plot(Xce,-FluxEz_x(index0,:,tindex),'r--');
hold on; plot(Xce,FluxEzce(index0,:,tindex),'g--','displayName','Flux_c_e');
lg2 = legend('show'); set(lg2,'location','northwest');
axis([-2*dX 3*dX 0 0.0002]);


%%%   compute Jz00 = d(rB)/dr/r
%
Jz00 = zeros(size(By));
for i=3:length(Xcc)-2
    Jz00(:,i,:) = (FluxEzce(:,i,:) - FluxEzce(:,i-1,:))/dX/Xcc(i);
end
Jz00(:,2,:) = Jz00(:,3,:);


figure(3); 
hold on; plot(Xcc,Jz00(index0,:,tindex),'Marker','*');
xlim([-2*dX 3*dX]); grid on;

