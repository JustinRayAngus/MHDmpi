%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look decomposition of U from 2D pinch module of m=0 mode
%%%   i.e. ... ExB, diamagnetic...polarization
%%%
%%%   E + UxB - Li[JxB - nabla(Pe)]/N = 0
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

numProcs = 2;
filePath = '../physicsMods/pinch2D/';
newDeck=0;
Li = 5.0e-4;

%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v4/test5_Li=4.6e-3_stable/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v4/test5_Li=4.6e-3_nuTherm/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v4/kR1/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/Hall_v4/test5_Li=5.0e-4_nuTherm_stable/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/kR10/Hall_v4/test5_Li=5.0e-4_stable/'; numProcs = 20; newDeck = 2;
%
filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_again/'; numProcs = 20; newDeck = 2;


thist = 6; %10.2; %9.4; %7.9;
thist2 = 6;  %10.2; %10.2; %9.4;

t0 = 3.6137e-8; % see normalizationParameters_zPinch.m
tA = 2.5e-8;    % from PoP 17, 072107 (2010)
%tA = 5.4036e-08;  % r0/VTi
%tA = 3.8209e-08;


tout = loadData(filePath,numProcs,'tout');
tout = tout*t0/tA; % normalize to 2010 paper alfven time
[~,tindex] = min(abs(tout-thist));
[~,tindex2] = min(abs(tout-thist2));

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

Jz0 = Jz;
Ez0 = Ez;
if(newDeck>=1)
  Jz0  = loadData(filePath,numProcs,'Jz0');
  Jx0  = loadData(filePath,numProcs,'Jx0');
  Ez0  = loadData(filePath,numProcs,'Ez0');
  Ex0  = loadData(filePath,numProcs,'Ex0');
  if(newDeck>=2)
      Li = loadData(filePath,numProcs,'lambda0');
      Se_source = loadData(filePath,numProcs,'Se_source');
      Si_source = loadData(filePath,numProcs,'Si_source');
      Pe  = loadData(filePath,numProcs,'Pe');
      Pi  = loadData(filePath,numProcs,'Pi');
      Te = loadData(filePath,numProcs,'Te');
      Ti = loadData(filePath,numProcs,'Ti');
  end
end


%%%   calculate pressure gradient
%
dPidx = zeros(size(Pi));
dPedx = zeros(size(Pe));
dPidz = zeros(size(Pi));
dPedz = zeros(size(Pe));
for i=2:length(Xcc)-1
    dPidx(:,i,:) = (Pi(:,i+1,:)-Pi(:,i-1,:))/(2*dX);
    dPedx(:,i,:) = (Pe(:,i+1,:)-Pe(:,i-1,:))/(2*dX);
end
for j=2:length(Zcc)-1
    dPidz(j,:,:) = (Pi(j+1,:,:)-Pi(j-1,:,:))/(2*dZ);
    dPedz(j,:,:) = (Pe(j+1,:,:)-Pe(j-1,:,:))/(2*dZ);
end
dPdz = dPidz + dPedz;
dPdx = dPidx + dPedx;


%%% calculate ideal electric field (E0 = -VxB)
%
Ex0 =  Vz.*By;
Ez0 = -Vx.*By;


%%% calculate non-ideal terms in Ohms law [E1 = Li*( JxB-nabla(Pe) )/N]
%
Ex1 = Li*(-Jz.*By - dPedx)./N;
Ez1 = Li*( Jx.*By - dPedz)./N;


%%%  calculate ideal ExB drift
%
V_ExB_z =  Ex./By;
V_ExB_x = -Ez./By;


%%% calculate non-ideal terms in velocity
%
Vz1 = Li*(Jz + dPedx./By)./N; 
Vx1 = Li*(Jx - dPedz./By)./N;


%%%   calculate diamagnetic drifts
%
%V_Di_z = Li*dPidx./By./N;
%V_Di_x = -Li*dPidz./By./N;


%%%   calculate polarization drifts
%
%V_Pi_z = Li*(Jz - dPdx./By)./N;
%V_Pi_x = Li*(Jx + dPdz./By)./N;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%           plot components of electric field
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f1=figure(1); 
set(f1,'position',[534 80 900 726]);

subplot(2,2,1);
dEx0 = Ex-Ex0;
pcolor(Zcc,Xcc,dEx0(:,:,tindex)'); 
shading flat; colorbar; box on;
xlabel('x direction'); ylabel('z direction');
title(['E_x-U_zB_y at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
subplot(2,2,3);
pcolor(Zcc,Xcc,Ex1(:,:,tindex)'); 
shading flat; colorbar; box on; 
caxis([min(min(dEx0(:,:,tindex))) max(max(dEx0(:,:,tindex)))]);
xlabel('x direction'); ylabel('z direction');
title(['(-J_zB_y - dP_e/dx)/N at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);


subplot(2,2,2);
dEz0 = Ez-Ez0;
pcolor(Zcc,Xcc,dEz0(:,:,tindex)'); 
shading flat; colorbar; box on;
xlabel('x direction'); ylabel('z direction');
title(['E_z+U_xB_y at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
subplot(2,2,4);
pcolor(Zcc,Xcc,Ez1(:,:,tindex)'); 
shading flat; colorbar; box on; 
caxis([min(min(dEz0(:,:,tindex))) max(max(dEz0(:,:,tindex)))]);
xlabel('x direction'); ylabel('z direction');
title(['(J_xB_y - dP_e/dz)/N at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%           plot components of ion velocity field
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f2=figure(2); 
set(f2,'position',[534 80 900 726]);

subplot(2,2,1);
dVx0 = Vx-V_ExB_x;
pcolor(Zcc,Xcc,dVx0(:,:,tindex)'); 
shading flat; colorbar; box on;
xlabel('x direction'); ylabel('z direction');
title(['U_x + E_z/B_y at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
subplot(2,2,3);
%Vx1 = V_Di_x + V_Pi_x;
pcolor(Zcc,Xcc,Vx1(:,:,tindex)'); 
shading flat; colorbar; box on; 
caxis([min(min(dVx0(:,:,tindex))) max(max(dVx0(:,:,tindex)))]);
xlabel('x direction'); ylabel('z direction');
title(['dP_i/dx/B_y/N at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);


subplot(2,2,2);
dVz0 = Vz-V_ExB_z;
pcolor(Zcc,Xcc,dVz0(:,:,tindex)'); 
shading flat; colorbar; box on;
xlabel('x direction'); ylabel('z direction');
title(['U_z - E_xB_y at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
subplot(2,2,4);
%Uz1 = U_Di_z + U_Pi_z;
pcolor(Zcc,Xcc,Vz1(:,:,tindex)'); 
shading flat; colorbar; box on; 
caxis([min(min(dVz0(:,:,tindex))) max(max(dVz0(:,:,tindex)))]);
xlabel('x direction'); ylabel('z direction');
title(['dP_i/dz/B_y/N at t = ',num2str(tout(tindex),3),' t_A']);
axis('square'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);


