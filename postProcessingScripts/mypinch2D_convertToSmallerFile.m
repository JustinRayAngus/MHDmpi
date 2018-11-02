%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  convert output variables from 2D pinch module of m=0 mode to
%%%  smaller file
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;



filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/ka10.0_nuTherm10/'; numProcs = 20; newDeck = 2;
% filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/ka3.0_nuTherm0/'; numProcs = 20; newDeck = 2;
% %filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/withShear/testing3/'; numProcs = 20; newDeck = 2;
% filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/withShear/ka3.0_M0.5_ideal/'; numProcs = 20; newDeck = 2;
% %filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/withShear/ka3.0_M0.5_taui1.0e-3_splitCspeed/'; numProcs = 20; newDeck = 2;
%

% filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/Paraschiv_Fig10/kR5.0_M0.85/'; numProcs = 20; newDeck = 2;
% filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/Paraschiv_Fig10/ideal_v0/kR5.0_M1.25/'; numProcs = 25; newDeck = 0;
% %filePath = '/Users/angus1/Documents/zPinch/myMHDsims/entropy_v1/Paraschiv_Fig10/ideal_v0/kR5.0_M2.5/'; numProcs = 50; newDeck = 0;
% 
% 
filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/Paraschiv_Fig10/kR5.0_M1.25/'; numProcs = 25; newDeck = 0;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/Paraschiv_Fig10/kR5.0_M1.25_320x800/'; numProcs = 100; newDeck = 0;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/Paraschiv_Fig10/kR5.0_M1.25_320x1200/'; numProcs = 120; newDeck = 0;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/Paraschiv_Fig10/kR5.0_M1.25_320x1600/'; numProcs = 144; newDeck = 0;

%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/Paraschiv_Fig10/kR5.0_M0.0_largerR/'; numProcs = 54; newDeck = 0;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/Paraschiv_Fig10/kR5.0_M0.0/'; numProcs = 25; newDeck = 0;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/Paraschiv_Fig10/kR5.0_M0.0_U1/'; numProcs = 25; newDeck = 0;

%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/minModBuild/testkR5/'; numProcs = 25; newDeck = 0;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/superBeeBuild/testkR5/'; numProcs = 25; newDeck = 0;
%filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/Paraschiv_Fig10/kR5.0_M0.0/'; numProcs = 25; newDeck = 0;

% filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/ka3/M0.5/'; numProcs = 25; newDeck = 0;
% filePath = '/Users/angus1/Documents/zPinch/myMHDsims/IdealMHD/ka3/M0.5_160x800/'; numProcs = 50; newDeck = 0;



tout = loadData(filePath,numProcs,'tout');
Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Zcc = loadData(filePath,numProcs,'Zcc');
Zce = loadData(filePath,numProcs,'Zce');
rcc = loadData(filePath,numProcs,'rcc');

N  = loadData(filePath,numProcs,'N');
S  = loadData(filePath,numProcs,'S');
By = loadData(filePath,numProcs,'By');
P  = loadData(filePath,numProcs,'P');
Ez  = loadData(filePath,numProcs,'Ez');
Ex  = loadData(filePath,numProcs,'Ex');
gamma0 = loadData(filePath,numProcs,'gamma0');



data.tout = tout;
data.Xcc = Xcc;
%data.Xce = Xce;
data.Zcc = Zcc;
%data.Zce = Zce;
data.rcc = rcc;
%
data.N = N;
data.S = S;
data.By = By;
data.Ez = Ez;
data.Ex = Ex;
%
data.gamma0 = gamma0;

save([filePath,'data.mat'],'data','-v7.3');







