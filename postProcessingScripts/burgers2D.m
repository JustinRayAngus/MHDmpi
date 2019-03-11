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
filePath = '../physicsMods/burgers2D/dataSave/';
filePath = '../physicsMods/burgers1D/';
%filePath = '../physicsMods/sodShock/';

Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Zcc = loadData(filePath,numProcs,'Zcc');
Zce = loadData(filePath,numProcs,'Zce');
tout= loadData(filePath,numProcs,'tout');
%
F0 = loadData(filePath,numProcs,'F0');
FluxLimL = loadData(filePath,numProcs,'FluxLimL_x');
FluxLimR = loadData(filePath,numProcs,'FluxLimR_x');
FluxL    = loadData(filePath,numProcs,'FluxL_x');
FluxR    = loadData(filePath,numProcs,'FluxR_x');
Flux     = loadData(filePath,numProcs,'Flux_x');


f1=figure(11); 
%set(f1,'position',[1030 925 1100 420]);
set(f1,'position',[1100 360 500 760]);
subplot(2,1,1);
hold on; pcolor(Xce(2:end-1),Zce(2:end-1),F0(3:end-1,3:end-1,1)); shading flat; colorbar;
xlabel('x direction'); caxis([0,1]);
ylabel('z direction'); box on;
title('initial profile');
axis('equal'); axis([-0.5 0.5 -0.5 0.5]);
%
subplot(2,1,2);
hold on; pcolor(Xce(2:end-1),Zce(2:end-1),F0(3:end-1,3:end-1,end)); colorbar;
xlabel('x direction');  caxis([0,1]);
ylabel('z direction'); shading flat; box on;
title('profile at final time step');
axis('equal'); axis([-0.5 0.5 -0.5 0.5]);
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
set(gca,'xtick',Xce(2):0.1:Xce(end-1));
xlabel('x'); ylabel('f');
title('function');
xlim([Xce(2) Xce(end-1)]);
%
subplot(2,1,2);
%hold on; plot(Xce,FluxL(:,1)); box on;
hold on; plot(Xce,FluxL(end/2,:,round(end/2)),'displayName','Left');
hold on; plot(Xce,FluxR(end/2,:,round(end/2)),'displayName','right');
%hold on; plot(Xce,FluxL(:,round(end))); grid on;
set(gca,'xtick',Xce(2):0.1:Xce(end-1));
xlabel('x'); ylabel('flux'); grid on;
title('flux'); xlim([Xce(2) Xce(end-1)]);
lg1=legend('show');

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


