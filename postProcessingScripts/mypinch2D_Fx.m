%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 2D pinch module of m=0 mode and 
%%%   investigate smoothing initial radial profile to achieve numerical
%%%   force balance
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/testingInitialFx/'; numProcs = 25;
filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/testingInitialFx2/'; numProcs = 100;

filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD_new/ka3.0/M0.35/80x400/'; numProcs = 50;

tout = loadData(filePath,numProcs,'tout');
%
Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Zcc = loadData(filePath,numProcs,'Zcc');
Zce = loadData(filePath,numProcs,'Zce');
rcc = loadData(filePath,numProcs,'rcc');
%
N  = loadData(filePath,numProcs,'N');
deltaP  = loadData(filePath,numProcs,'deltaP');
%P0  = hdf5read(thisFile,'P0');
Mx  = loadData(filePath,numProcs,'Mx');
Mz  = loadData(filePath,numProcs,'Mz');
S  = loadData(filePath,numProcs,'S');
By = loadData(filePath,numProcs,'By');
P  = loadData(filePath,numProcs,'P');
P0  = loadData(filePath,numProcs,'P0');
%
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


%%%   define Bennett profiles for P and B
%
By0 = sqrt(2)*Xcc./(1+Xcc.^2);
P0  = 1./(1+Xcc.^2).^2;


%%%   plot contours
%
close(figure(13));
f13=figure(13); set(f13,'position',[1120 300 1200 800]);
%
subplot(2,3,1);
plot(Xcc,N(3,:,1)); box on;
xlabel('x direction'); 
title('density');
axis('square'); axis([Xce(2) Xce(end-1) 0 1.2]);
%
subplot(2,3,2);
plot(Xcc,P(3,:,1)); box on;
hold on; plot(Xcc,P0,'r--');
xlabel('x direction'); 
title('pressure');
axis('square'); axis([Xce(2) Xce(end-1) 0 1.2]);
%
subplot(2,3,3);
plot(Xcc,By(3,:,1)); box on;
hold on; plot(Xcc,By0,'r--');
xlabel('x direction'); 
title('magnetic field');
axis('square'); axis([Xce(2) Xce(end-1) 0 1.2]);
%
subplot(2,3,4);
plot(Xcc,Fx(3,:,1)); box on;
xlabel('x direction'); 
title('radial force');
axis('square'); xlim([Xce(2) Xce(end-1)]);
%
subplot(2,3,5);
plot(Xcc,P(3,:,1)-P0'); box on;
xlabel('x direction'); 
title('P-P0');
axis('square'); xlim([Xce(2) Xce(end-1)]);
%
subplot(2,3,6);
plot(Xcc,By(3,:,1)-By0'); box on;
xlabel('x direction'); 
title('By-By0');
axis('square'); xlim([Xce(2) Xce(end-1)]);

