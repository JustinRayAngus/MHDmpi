%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from burgers2D simulation
%%%
%%%   superbee works better than Van leer
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

numProcs = 4;
filePath = '../physicsMods/burgers2D/';
%filePath = '../physicsMods/sodShock/';

for i=1:numProcs
fileName = ['output',num2str(i-1),'.h5'];
thisFile = [filePath,fileName];
procID  = hdf5read(thisFile,'procID');
fileinfo = hdf5info(thisFile);
Xcc = hdf5read(thisFile,'Xcc');
Zcc = hdf5read(thisFile,'Zcc');
Xce = hdf5read(thisFile,'Xce');


F0  = hdf5read(thisFile,'F0');
FluxLimL  = hdf5read(thisFile,'FluxLimL_x');
FluxLimR  = hdf5read(thisFile,'FluxLimR_x');
FluxL    = hdf5read(thisFile,'FluxL_x');
FluxR    = hdf5read(thisFile,'FluxR_x');
Flux  = hdf5read(thisFile,'Flux_x');
tout= hdf5read(thisFile,'tout');


f11=figure(12); 
%set(f1,'position',[1030 925 1100 420]);
set(f11,'position',[1800 360 500 760]);
subplot(2,1,1);
hold on; pcolor(Xcc,Zcc,F0(:,:,1)); shading flat;colorbar;
xlabel('x direction'); axis('equal'); caxis([0,1]);
ylabel('z direction');
title('initial profile');
%
subplot(2,1,2);
hold on; pcolor(Xcc,Zcc,F0(:,:,round(end))); colorbar;
xlabel('x direction'); axis('equal'); caxis([0,1]);
ylabel('z direction'); shading flat;
title('profile at final time step');
%
map = colormap('jet');
map(1,:) = [1 1 1];
colormap(map)
    
%%%
%
%close(figure(1));
f1=figure(2); 
%set(f1,'position',[1030 925 1100 420]);
set(f1,'position',[1190 360 500 760]);


subplot(2,1,1);
hold on; plot(Xcc,F0(end/2,:,1),'black'); box on;
hold on; plot(Xcc,F0(end/2,:,round(end/2)),'b');
hold on; plot(Xcc,F0(end/2,:,end),'r'); grid on;
%hold on; plot(Zcc,F0(:,end/2,1)); box on;
%hold on; plot(Zcc,F0(:,end/2,round(end/2)));
%hold on; plot(Zcc,F0(:,end/2,end)); grid on;
set(gca,'xtick',-0.5:0.1:0.5);
xlabel('x'); ylabel('f');
title('function');
xlim([-0.5 0.5]);
%
subplot(2,1,2);
%hold on; plot(Xce,FluxL(:,1)); box on;
hold on; plot(Xce,FluxL(end/2,:,round(end/2)),'displayName','Left');
hold on; plot(Xce,FluxR(end/2,:,round(end/2)),'displayName','right');
%hold on; plot(Xce,FluxL(:,round(end))); grid on;
set(gca,'xtick',-0.5:0.1:0.5);
xlabel('x'); ylabel('flux'); grid on;
title('flux'); xlim([-0.5 0.5]);
if(i==1) 
    lg1=legend('show');
end

%
%
%

% f3=figure(3); set(f3,'position',[1030 430 1100 420]);
% %
% subplot(1,2,1);
% hold on; plot(Xce,FluxRatio(:,end)); grid on;
% set(gca,'xtick',-0.5:0.1:0.5); box on;
% xlim([-0.5 0.5]);
% %
% subplot(1,2,2);
% hold on; plot(Xce,FluxLim(:,end)); grid on; box on;
% set(gca,'xtick',-0.5:0.1:0.5);
% xlim([-0.5 0.5]);

%
%
%

% f4=figure(4); set(f4,'position',[1030 0 1100 420]);
% %
% subplot(1,2,1);
% hold on; plot(Xce,FluxL(:,end)); grid on;
% set(gca,'xtick',-0.5:0.1:0.5); box on;
% title('left going Flux');
% xlim([-0.5 0.5]);
% %
% subplot(1,2,2);
% hold on; plot(Xce,FluxR(:,end)); grid on; box on;
% set(gca,'xtick',-0.5:0.1:0.5);
% title('right going Flux');
% xlim([-0.5 0.5]);

end

