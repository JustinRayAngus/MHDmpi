%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  make movie from 2D pinch module of m=0 mode
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
%filePath = '../../fromQuartz/pinch2D/kR10/ideal2/testing5_stable/'; numProcs = 20;


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
filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/test5_Li=0.0/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/test5_Li=5.0e-3/'; numProcs = 20; newDeck = 1;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v2/test5_Li=1.0e-1/'; numProcs = 20; newDeck = 1;
%
filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_again/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_tauei/'; numProcs = 20; newDeck = 2;
% filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable_tauei/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-2_again/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=2.0e-2_again/'; numProcs = 20; newDeck = 2;

filePath = '../../fromQuartz/pinch2D/entropy_v0/kR10/test5_Li=5.0e-3_stable/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v0/kR10/testing_Li0.0/'; numProcs = 20; newDeck = 2;

filePath = '../../fromQuartz/pinch2D/entropy_v0/Li1.5e-2/ka3.0/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v0/Li5.0e-2/ka20_smallerDt/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v0/Li0.12/ka10.0/'; numProcs = 20; newDeck = 2;

%filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/ka3.0/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/ka3.0_noGyroVisc/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/noGyroVisc/ka10.0_taui1.0e-2/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/noGyroVisc/ka3.0_taui1.0e-2_nuTherm100/'; numProcs = 20; newDeck = 2;


%t0 = 3.6137e-8; % see normalizationParameters_zPinch.m
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
Jx  = loadData(filePath,numProcs,'Jx');
Cs  = loadData(filePath,numProcs,'Cs');
%Te  = loadData(filePath,numProcs,'Te');
gamma0 = loadData(filePath,numProcs,'gamma0');
Li0 = loadData(filePath,numProcs,'lambda0');
dX = Xcc(2)-Xcc(1);
dZ = Zcc(2)-Zcc(1);
Fx = loadData(filePath,numProcs,'Fx'); % - d(P+B^2/2)/dr - B^2/r

T = P./N;
Ptot = P+By.^2/2+Mx.*Vx/2.0;


%%%   plot contours
%
close(figure(1));
f1=figure(1); set(f1,'position',[1100 300 1500 500]);
%f1=figure(1); set(f1,'position',[400 400 1800 200]);
set(gcf,'color','w');

thisFigNum = 0;
for thist=1:length(tout)-1

    subplot(1,3,1);
    h1=pcolor(Zcc,Xcc,N(:,:,thist)'); colorbar; box on
    xlabel('z/r_0'); shading flat;
    ylabel('r/r_0'); 
    title(['density at t = ',num2str(tout(thist),3),' t_0=c_s/r_0']);
    axis('equal'); 
    axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
    %
    %
    %
    subplot(1,3,2);
    h2=pcolor(Zcc,Xcc,By(:,:,thist)'); colorbar; box on
    xlabel('z/r_0'); shading flat;
    ylabel('r/r_0'); 
    title(['B_y at t = ',num2str(tout(thist),3),' t_0=c_s/r_0']);
    axis('equal'); 
    axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
    colormap('jet');
    %
    %
    %
    subplot(1,3,3);
   % h3=pcolor(Zcc,Xcc,Te(:,:,thist)'); colorbar; box on
    h3=pcolor(Zcc,Xcc,Ez(:,:,thist)'); colorbar; box on
    xlabel('z/r_0'); shading flat;
    ylabel('r/r_0'); 
    title(['E_z at t = ',num2str(tout(thist),3),' t_0=c_s/r_0']);
    axis('equal'); 
    axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
    colormap('jet');    
    
%     %%%   put time stamp on figure
%     %
%     a1=annotation(f1,'textbox',...
%     [0.03 0.5 0.1 0.05],...
%     'String',['t = ',num2str(tout(thist)),' t_A'],...
%     'FitBoxToText','off', ...
%     'backgroundcolor','y');
    
    
    M(thisFigNum+1) = getframe(f1); 
    thisFigNum = thisFigNum+1;
    if(tout(thist)~=tout(end-1))
      %  delete(a1);
        delete(h1);
        delete(h2);
        delete(h3);
    end

end

f1pos = get(f1,'position');
f11=figure(11);
set(f11,'Position',f1pos);
movie(f11,M,1,3);
%v=VideoWriter([dataPath,'movieFigs/2DrBthMovie.avi']);
v=VideoWriter('./2Dpinch.avi');
v.FrameRate = 4; %v.Quality=100; %v.CompressionRatio = 2;
%v.LosslessCompression = true;
open(v); writeVideo(v,M);
close(v);
