%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from MHD1D simulation
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

numProcs = 2;
filePath = '../build/';
%filePath = '../../MHD1D_MPI_save/buildMPI/';

for i=1:numProcs
fileName = ['output',num2str(i-1),'.h5'];
thisFile = [filePath,fileName];
fileinfo = hdf5info(thisFile);
Xcc = hdf5read(thisFile,'Xcc');
Xce = hdf5read(thisFile,'Xce');
F0  = hdf5read(thisFile,'F0');
Flux  = hdf5read(thisFile,'Flux');
tout= hdf5read(thisFile,'tout');


%%%
%
%close(figure(1));
f1=figure(1); set(f1,'position',[646     1  1100   420]);
%
subplot(1,2,1);
hold on; plot(Xcc,F0(:,1)); box on;
hold on; plot(Xcc,F0(:,round(end/2)));
hold on; plot(Xcc,F0(:,end)); grid on;
xlim([0 1]);
%
subplot(1,2,2);
hold on; plot(Xce,Flux(:,1)); box on;
hold on; plot(Xce,Flux(:,round(end/2)));
hold on; plot(Xce,Flux(:,round(end))); grid on;
xlim([0 1]);
end


%figure(3); hold on; plot(Xcc); box on;