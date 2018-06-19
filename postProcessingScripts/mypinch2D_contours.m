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
newDeck=1;


%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/testing7/'; numProcs = 20; newDeck = 1;
filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/test5_Li=1.0e-1/'; numProcs = 20; newDeck = 1;
lambdai = 1.0e-1;
thist = 3.6; %10.2; %9.4;

% filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/test5_Li=0.0/'; numProcs = 20; newDeck = 1;
% lambdai = 4.6e-3;
% thist = 7.0; %9.4; %7.9;


t0 = 3.6137e-8; % see normalizationParameters_zPinch.m
tA = 2.5e-8;    % from PoP 17, 072107 (2010)
%tA = 5.4036e-08;  % r0/VTi
%tA = 3.8209e-08;


tout = loadData(filePath,numProcs,'tout');
tout = tout*t0/tA; % normalize to 2010 paper alfven time
[~,tindex] = min(abs(tout-thist));
%[~,tindex2] = min(abs(tout-thist2));

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

FluxEz_x  = loadData(filePath,numProcs,'FluxEz_x');
FluxEz_x2 = loadData(filePath,numProcs,'FluxEz_x2');
FluxEx_z  = loadData(filePath,numProcs,'FluxEx_z');
FluxEx_z2 = loadData(filePath,numProcs,'FluxEx_z2');

Jx0  = loadData(filePath,numProcs,'Jx0');
Jx02  = loadData(filePath,numProcs,'Jx02');

Vx  = loadData(filePath,numProcs,'Vx');
Vz  = loadData(filePath,numProcs,'Vz');
Jz  = loadData(filePath,numProcs,'Jz');
Jx  = loadData(filePath,numProcs,'Jx');
Ez  = loadData(filePath,numProcs,'Ez');
Ez0  = loadData(filePath,numProcs,'Ez0');
Ezh = lambdai*Jx.*By./N;
Ex  = loadData(filePath,numProcs,'Ex');
Cs  = loadData(filePath,numProcs,'Cs');
gamma0 = loadData(filePath,numProcs,'gamma0');
dX = Xcc(2)-Xcc(1);
dZ = Zcc(2)-Zcc(1);


%%%   plot contours
%
close(figure(1));
f1=figure(1); set(f1,'position',[1637 50 1200 1244]);
%
subplot(3,2,1);
pcolor(Zcc,Xcc,N(:,:,tindex)'); colorbar; box on
xlabel('z direction'); shading flat;
ylabel('r direction');
title(['density at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
subplot(3,2,2);
pcolor(Zcc,Xcc,By(:,:,tindex)'); colorbar; box on
xlabel('z direction'); shading flat;
ylabel('r direction');
title(['By at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
FluxEx_diff = abs(FluxEx_z-FluxEx_z2);
subplot(3,2,4);
pcolor(Zcc,Xcc,Ezh(:,:,tindex)'); colorbar; box on
xlabel('z direction'); shading flat;
ylabel('r direction');
title(['hall Ez at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
subplot(3,2,3);
pcolor(Zcc,Xcc,Ez0(:,:,tindex)'); colorbar; box on
xlabel('z direction'); shading flat;
ylabel('r direction');
title(['ideal Ez at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
subplot(3,2,5);
pcolor(Zcc,Xcc,Ez(:,:,tindex)'); colorbar; box on
xlabel('z direction'); shading flat;
ylabel('r direction');
title(['total Ez at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
FluxEz_diff = abs(FluxEz_x-FluxEz_x2);
subplot(3,2,6);
pcolor(Zcc,Xcc,Jx(:,:,tindex)'); colorbar; box on
xlabel('z direction'); shading flat;
ylabel('r direction');
title(['Jx at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);


map = colormap('jet');
%map(1,:) = [1 1 1];
colormap(map)

%
%
%



