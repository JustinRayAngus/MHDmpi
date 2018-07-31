%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  make figures for ICOPS poster
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
%filePath = '../../fromQuartz/pinch2D/data_200nx_kR1/';
%filePath = '../../fromQuartz/pinch2D/data_200nx_kR5/';
filePath = '../../fromQuartz/pinch2D/data_200nx_kR10/';
filePath = '../../fromQuartz/pinch2D/data_400nx_kR10/'; numProcs = 20;
filePath = '../../fromQuartz/pinch2D/data_800nx_kR10/'; numProcs = 40;
filePath = '../../fromQuartz/pinch2D/data_160nz_800nx_kR10/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/data_200nx_kR10_ideal1/';


filePath = '../../fromQuartz/pinch2D/kR1/data_160nz_800nx/'; numProcs = 80;
filePath = '../../fromQuartz/pinch2D/kR1/data_320nz_800nx/'; numProcs = 80;

%filePath = '../../fromQuartz/pinch2D/kR5/data_160nz_800nx/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR5/data_320nz_800nx/'; numProcs = 80;

filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing1/'; numProcs = 20;
filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing4/'; numProcs = 20;
filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing5/'; numProcs = 20;
filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing5_stable/'; numProcs = 20;


% filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test5_Li=0/'; numProcs = 20; newDeck = 1;
% filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test5_Li=4.6e-3/'; numProcs = 20; newDeck = 1;
% filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test5_Li=1.0e-2/'; numProcs = 20; newDeck = 1;
% %filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test4_Li=1.0e-2/'; numProcs = 20; newDeck = 1;
% 
% %filePath = '../../fromQuartz/pinch2D/kR10/Hall_v0/test6_Li=1.0e-3/'; numProcs = 20; newDeck = 1;
% %filePath = '../../fromQuartz/pinch2D/kR10/Hall_v1/test4_Li=1.0e-3/'; numProcs = 20; newDeck = 1;
% %filePath = '../../fromQuartz/pinch2D/kR10/Hall_v1/test4_Li=1.0e-2/'; numProcs = 20; newDeck = 1;
% filePath = '../../fromQuartz/pinch2D/kR10/Hall_v1/test4_Li=1.0e-1/'; numProcs = 20; newDeck = 1;
% %filePath = '../../fromQuartz/pinch2D/kR10/Hall_v1/test5_Li=1.0e-2/'; numProcs = 20; newDeck = 1;
% %
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/test5_Li=0.0/'; numProcs = 20; newDeck = 1;
% %filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/test5_Li=1.0e-1/'; numProcs = 20; newDeck = 1;
% filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/test5_Li=4.6e-3/'; numProcs = 20; newDeck = 1;
% filePath = '../../fromQuartz/pinch2D/kR10/Hall_v3/test5_Li=4.6e-3_stable/'; numProcs = 20; newDeck = 1;
%
% filePath = '../../fromQuartz/pinch2D/kR10/Hall_v3/test5_Li=1.0e-4_stable/'; numProcs = 20; newDeck = 2;
%
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_tauei/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable_tauei/'; numProcs = 20; newDeck = 2;
%%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-2/'; numProcs = 20; newDeck = 2;


filePath = '../../fromQuartz/pinch2D/entropy_v0/Li0.0/ka3.0/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v0/Li1.5e-2/ka3.0_stable_nuT/'; numProcs = 20; newDeck = 2;



t0 = 1.2046e-8; % using pinch radius for length scale
tA = 2.5e-8;    % from PoP 17, 072107 (2010)


tout = loadData(filePath,numProcs,'tout');
%tout = tout*t0/tA; % normalize to 2010 paper alfven time
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
%Te  = loadData(filePath,numProcs,'Te');
gamma0 = loadData(filePath,numProcs,'gamma0');
dX = Xcc(2)-Xcc(1);
dZ = Zcc(2)-Zcc(1);
Fx = loadData(filePath,numProcs,'Fx'); % - d(P+B^2/2)/dr - B^2/r

T = P./N;
Ptot = P+By.^2/2+Mx.*Vx/2.0;

%tout0 = [7.9 9.4 10.2]*tA/t0; % output times in unites of t0/tA
tout0 = [6 25 37.6];      % for drift-ideal with gamma=2.1
[~,itout0] = min(abs(tout-tout0));
if(itout0(end)==length(tout))
    itout0(end) = itout0(end)-1;
end


%%%   plot contours
%
close(figure(2));
f1=figure(2); set(f1,'position',[1300 80 1000 1220]);
%f1=figure(1); set(f1,'position',[400 400 1800 200]);
set(gcf,'color','w');


subplot(3,3,1);
h1=pcolor(Zcc,Xcc,N(:,:,itout0(1))'); colorbar; box on
xlabel('z/a'); shading flat; caxis([0 1.4]);
ylabel('r/a'); 
title(['N at t = ',num2str(tout(itout0(1)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(3,3,2);
h1=pcolor(Zcc,Xcc,N(:,:,itout0(2))'); colorbar; box on
xlabel('z/a'); shading flat; caxis([0 1.4]);
ylabel('r/a'); 
title(['N at t = ',num2str(tout(itout0(2)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(3,3,3);
h1=pcolor(Zcc,Xcc,N(:,:,itout0(3))'); colorbar; box on
xlabel('z/a'); shading flat; caxis([0 1.4]);
ylabel('r/a'); 
title(['N at t = ',num2str(tout(itout0(3)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
subplot(3,3,4);
h2=pcolor(Zcc,Xcc,By(:,:,itout0(1))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['B_y at t = ',num2str(tout(itout0(1)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(3,3,5);
h2=pcolor(Zcc,Xcc,By(:,:,itout0(2))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['B_y at t = ',num2str(tout(itout0(2)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(3,3,6);
h2=pcolor(Zcc,Xcc,By(:,:,itout0(3))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['B_y at t = ',num2str(tout(itout0(3)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
subplot(3,3,7);
% h3=pcolor(Zcc,Xcc,Te(:,:,thist)'); colorbar; box on
h3=pcolor(Zcc,Xcc,Ez(:,:,itout0(1))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['E_z at t = ',num2str(tout(itout0(1)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(3,3,8);
% h3=pcolor(Zcc,Xcc,Te(:,:,thist)'); colorbar; box on
h3=pcolor(Zcc,Xcc,Ez(:,:,itout0(2))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['E_z at t = ',num2str(tout(itout0(2)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(3,3,9);
% h3=pcolor(Zcc,Xcc,Te(:,:,thist)'); colorbar; box on
h3=pcolor(Zcc,Xcc,Ez(:,:,itout0(3))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['E_z at t = ',num2str(tout(itout0(3)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
colormap('jet');    
%colormap('hot'); 
    


