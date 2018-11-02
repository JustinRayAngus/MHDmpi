%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 2D pinch module of m=0 mode and 
%%%   investigate source of asymmetry for certain resolutions
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




%%%   kR = 10
%
filePath = '../../fromQuartz/pinch2D/kR10/data_80nz_400nx/'; numProcs = 20;
filePath = '../../fromQuartz/pinch2D/kR10/data_160nz_800nx/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/data_400nz_800nx/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/testing/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/kR10/testing2/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/testing3/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/testing4/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR10/data_80nz_400nx_ideal1/'; numProcs = 20;
%

filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/ka3/M0.0/160x800_lasttry/'; numProcs = 100;

plotBackIndex = 1; % plot time will be end-plotBackIndex
thist = 1; %9.4; %7.9;
thist2 = 2; %10.2; %9.4;

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
T = P./N;
Te = T/2;
Ti = T/2;
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



%%%   compute z-symmetric error for density
%
errorN = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
errorS = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
errorP = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
errorB = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
for j=1:length(errorN(:,1,1))
    errorN(j,:,:) = N(j+2,:,:)  - N(end-1-j,:,:);
    errorB(j,:,:) = By(j+2,:,:) - By(end-1-j,:,:);
    errorS(j,:,:) = S(j+2,:,:)  - S(end-1-j,:,:);
    errorP(j,:,:) = P(j+2,:,:)  - P(end-1-j,:,:);
end

figure(7); 
pcolor(Zcc(3:length(errorN(:,1,1))+2),Xcc,log10(abs(errorN(:,:,2)))'); 
shading flat; colorbar; title('error in N symmetry');

figure(8); 
pcolor(Zcc(3:length(errorN(:,1,1))+2),Xcc,log10(abs(errorP(:,:,2)))'); 
shading flat; colorbar; title('error in P symmetry');

colormap('jet');

maxErrorN0 = max(max(errorN(:,:,1)))
maxErrorP0 = max(max(errorP(:,:,1)))

