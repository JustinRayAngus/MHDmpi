%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  look at movies for electro-thermal instability sims
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
mu0 = 4*pi*1e-7;
cvac = 2.9979e10; % speed of light [cm/s]
me = 9.1094e-28;  % electron mass [g]
Mi = 2.0*1.6605e-24;    % ion mass [g]
meMi = me/Mi;


numProcs = 4;
filePath = '../physicsMods/liftOff2Drth/data0/';
filePath = '../physicsMods/liftOff2Drth/';

%procID  = hdf5read(thisFile,'procID');
%fileinfo = hdf5info(thisFile);
%
Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Ycc = loadData(filePath,numProcs,'Zcc');
Yce = loadData(filePath,numProcs,'Zce');

Bx  = loadData(filePath,numProcs,'Bx');
By  = loadData(filePath,numProcs,'By');
Bz  = loadData(filePath,numProcs,'Bz');
%
Ex  = loadData(filePath,numProcs,'Ex');
Ey  = loadData(filePath,numProcs,'Ey');
Ez  = loadData(filePath,numProcs,'Ez');
%
Ez_cc  = loadData(filePath,numProcs,'Ez_cc');
%
Ee  = loadData(filePath,numProcs,'Ee');
Te  = loadData(filePath,numProcs,'Te');
taue  = loadData(filePath,numProcs,'taue');
SourceEe  = loadData(filePath,numProcs,'SourceEe');
etaJsq = loadData(filePath,numProcs,'etaJsq');
Stherm = loadData(filePath,numProcs,'Stherm');

qex  = loadData(filePath,numProcs,'qex');
qey  = loadData(filePath,numProcs,'qey');
qex0  = loadData(filePath,numProcs,'qex0');
qey0  = loadData(filePath,numProcs,'qey0');

sig  = loadData(filePath,numProcs,'sig');
sigce_x  = loadData(filePath,numProcs,'sigce_x');
sigce_y  = loadData(filePath,numProcs,'sigce_y');
sigce_xy  = loadData(filePath,numProcs,'sigce_xy');
Jz = 1e9*sigce_xy.*Ez; % [statA/cm^2]
Jz_kApcmsq = Jz/(3e12); % [kA/cm^2]

%
gamma0 = loadData(filePath,numProcs,'gamma0');
N0 = loadData(filePath,numProcs,'N0');
Ti0 = loadData(filePath,numProcs,'Ti0');
%
tns = loadData(filePath,numProcs,'tout'); % [ns]
%
dX = Xcc(2)-Xcc(1);
dY = Ycc(2)-Ycc(1); 

%%% calculate Jz=cvac/(4*pi)*curlBy
%
curlBy = zeros(length(Xce),length(tns));
for i=2:length(Xcc)-2
    curlBy(i,:) = (Xcc(i+1)*By(16,i+1,:) ...
                -  Xcc(i)*By(16,i,:))/dX/Xce(i);
end
Jz0 = curlBy/4/pi*cvac; % [statA/cm^2]
%Jz0(end-1,:) = 0;
Jz0_kApcmsq = Jz0/(3e12); % [kA/cm^2]

Stherm2 = 2*meMi./taue/(gamma0-1)*N0.*(Te-Ti0);
deltaT_HeatBalance = etaJsq/(2.0*meMi*N0)*(gamma0-1).*taue;
Te_HeatBalance = deltaT_HeatBalance + Ti0;

thisit = length(tns);

close(figure(1));
f1=figure(1); set(f1,'position',[1000 960 1600 400]);
%
subplot(1,3,1);
h11=pcolor(Ycc,Xce,Bx(:,:,thisit)'); shading flat; colorbar; box on;
xlabel('y [cm]'); ylabel('x [cm]'); 
title(['B_x at t = ',num2str(tns(thisit),3),' ns']);
axis('equal'); 
axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);
%
subplot(1,3,2);
h12=pcolor(Yce,Xcc,By(:,:,thisit)'); shading flat; colorbar; box on;
xlabel('y [cm]'); ylabel('x [cm]'); 
title(['B_y at t = ',num2str(tns(thisit),3),' ns']);
axis('equal'); 
axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);
%
subplot(1,3,3);
h12=pcolor(Ycc,Xcc,Bz(:,:,thisit)'); shading flat; colorbar; box on;
xlabel('y [cm]'); ylabel('x [cm]'); 
title(['B_z at t = ',num2str(tns(thisit),3),' ns']);
axis('equal'); 
axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);


close(figure(2));
f2=figure(2); set(f2,'position',[1000 360 1600 400]);
%
subplot(1,3,1);
h21=pcolor(Yce,Xcc,Ex(:,:,thisit)'); shading flat; colorbar; box on;
xlabel('y [cm]');
ylabel('x [cm]'); 
title(['E_x at t = ',num2str(tns(thisit),3),' ns']);
axis('equal'); 
axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);
%
subplot(1,3,2);
h22=pcolor(Ycc,Xce,Ey(:,:,thisit)'); shading flat; colorbar; box on;
xlabel('y [cm]');
ylabel('x [cm]'); 
title(['E_y at t = ',num2str(tns(thisit),3),' ns']);
axis('equal'); 
axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);
%
subplot(1,3,3);
h23=pcolor(Yce,Xce,Ez(:,:,thisit)'); shading flat; colorbar; box on;
xlabel('y [cm]');
ylabel('x [cm]'); 
title(['E_z at t = ',num2str(tns(thisit),3),' ns']);
axis('equal'); 
axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);


f3=figure(3);
h3=pcolor(Ycc,Xcc,Te(:,:,thisit)'); shading flat; colorbar; box on;
xlabel('y [cm]');
ylabel('x [cm]'); 
title(['T_e at t = ',num2str(tns(thisit),3),' ns']);
axis('equal'); 
axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);


%%%   compare Te with that from energy balance
%
f4=figure(4);
plot(tns,squeeze(Te(16,3,:)),'displayName','T_e');
box on; grid on;
xlabel('time [ns]');
ylabel('boundary Te [eV]');


%%%   compare Jz with Jz0
%
f5=figure(5);
plot(Xce,Jz_kApcmsq(16,:,end),'displayName','J_z');
hold on; plot(Xce,Jz0_kApcmsq(:,end),'r--', ...
              'displayName','from curlBy');
box on; grid on;
legend('show','location','best');

%%%   compare Jz with Jz0
%
f6=figure(6);
plot(Xcc,Te(16,:,end),'displayName','T_e');
hold on; plot(Xcc,Te_HeatBalance(16,:,end),'r--', ...
              'displayName','from heat balance');
box on; grid on;
legend('show','location','best');


%%%   plot maximum Bx to see what growth rate is
%
Bxmax = zeros(size(tns));
for i=1:length(tns)
    Bxmax(i) = max(max(Bx(:,:,i)));
end

y = log(Bxmax);
gamma = zeros(size(tns));
for i=1:length(tns)-1
    gamma(i) = (y(i+1)-y(i))/(tns(i+1)-tns(i));
end


f7=figure(7); 
%
subplot(1,2,1);
plot(tns,log(Bxmax));
title('log maximum B_x');
xlabel('time [ns]');
box on; grid on;

subplot(1,2,2);
plot(tns,gamma); box on; grid on;
title('growth rate [nHz]');
ylabel('time [ns]');
