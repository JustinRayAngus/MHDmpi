%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  look at electric field from m=0 simulations
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

%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_tauei/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-3_stable_tauei/'; numProcs = 20; newDeck = 2;
%%filePath = '../../fromQuartz/pinch2D/kR10/entropy_v0/test5_Li=5.0e-2/'; numProcs = 20; newDeck = 2;


filePath = '../../fromQuartz/pinch2D/entropy_v0/Li0.0/ka3.0/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/ka3.0_noGyroVisc_taui1.0e-2/'; numProcs = 20; newDeck = 2;



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
%
Jx0 = loadData(filePath,numProcs,'Jx0'); % curlB dot x
Jz0 = loadData(filePath,numProcs,'Jz0'); % curlB dot z
%
%Te  = loadData(filePath,numProcs,'Te');
gamma0 = loadData(filePath,numProcs,'gamma0');
dX = Xcc(2)-Xcc(1);
dZ = Zcc(2)-Zcc(1);
Fx = loadData(filePath,numProcs,'Fx'); % - d(P+B^2/2)/dr - B^2/r

T = P./N;
Ptot = P+By.^2/2+Mx.*Vx/2.0;

%tout0 = [7.9 9.4 10.2]*tA/t0; % output times in unites of t0/tA
tout0 = [6 14.8 20];      % for ideal MHD
tout0 = [6 10 14.8];    % for drift-ideal MHD
[~,itout0] = min(abs(tout-tout0));
if(itout0(end)==length(tout))
    itout0(end) = itout0(end)-1;
end


%%%   plot contours
%
close(figure(1));
f1=figure(1); set(f1,'position',[1300 80 1000 1220]);
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
    

%%%   calculate divergence and curl of electric field
%
divE = zeros(size(N)); % charge density
dEzdz_cez = zeros(length(Zce),length(Xcc),length(tout));
dExdz_cez = zeros(length(Zce),length(Xcc),length(tout));
dEzdx_cez = zeros(length(Zce),length(Xcc),length(tout));
drErdrr_cez = zeros(length(Zce),length(Xcc),length(tout));
drErdrr_cc = zeros(length(Zcc),length(Xce),length(tout));
%divE_ce = zeros(length(Zce),length(Xce),length(tout));
curlE = zeros(size(N)); % only in y-direction
for i=2:length(Xcc)-1
    for j=2:length(Zcc)-1
        divE(j,i,:) = (Ez(j+1,i,:)-Ez(j-1,i,:))/2/dZ ...
                    + (Xcc(i+1)*Ex(j,i+1,:)-Xcc(i-1)*Ex(j,i-1,:))/2/dX/Xcc(i);
        %
        dEzdz_cez(j,i,:)  = (Ez(j+1,i,:)-Ez(j,i,:))/dZ;
        dExdz_cez(j,i,:)  = (Ex(j+1,i,:)-Ex(j,i,:))/dZ;
        drErdrr_cc(j,i,:) = (Xcc(i+1)*Ex(j,i+1,:)-Xcc(i-1)*Ex(j,i-1,:))/2/dX/Xcc(i);
        drErdrr_cez(j,i,:) = ( Xcc(i+1)*(Ex(j,i+1,:) + Ex(j+1,i+1,:))/2 ...
                           -   Xcc(i-1)*(Ex(j,i-1,:) + Ex(j+1,i-1,:))/2 )/2/dX/Xcc(i);
        dEzdx_cez(j,i,:) = ( (Ez(j,i+1,:) + Ez(j+1,i+1,:))/2 ...
                           - (Ez(j,i-1,:) + Ez(j+1,i-1,:))/2 )/2/dX;
        %
        curlE(j,i,:) = (Ex(j+1,i,:)-Ex(j-1,i,:))/2/dZ ...
                     - (Ez(j,i+1,:)-Ez(j,i-1,:))/2/dX;
    end
end
divE_cez = dEzdz_cez + drErdrr_cez;
divE_cez(1,:,:)    = divE_cez(end-2,:,:);
curlE_cez = dExdz_cez - dEzdx_cez;
curlE_cez(1,:,:)    = curlE_cez(end-2,:,:);
%divE_cez(end,:,:)  = divE(4,:,:);
divE(1,:,:)    = divE(end-3,:,:);
divE(end,:,:)  = divE(4,:,:);
curlE(1,:,:)   = curlE(end-3,:,:);
curlE(end,:,:) = curlE(4,:,:);


nX = length(Xcc);
nZ = length(Zcc);
xcc = zeros(nX,nZ);
zcc = zeros(nX,nZ);
for j=1:nZ
    xcc(:,j) = Xcc';
end
for i=1:nX
    zcc(i,:) = Zcc';
end


phi = zeros(size(N)); % electrostatic potential
Cy = zeros(size(N));  % vector potential of E field (See Helmotlz Decomp)
Ax  = zeros(size(N)); % electromagnetic potential
Az  = zeros(size(N)); % electromagnetic potential
usePoisson = 1;
if(usePoisson)
    
    %%%   calculate potentials by solving poissons equation
    %
    for n=1:length(tout)
        phi0 = Poisson_Solver_2D(xcc,zcc,-divE(:,:,n)',0,1);
        phi(:,:,n) = phi0';
        %
        Ax0 = Poisson_Solver_2D(xcc,zcc,-Jx0(:,:,n)',0,1);
        Ax(:,:,n) = Ax0';
        %
        Az0 = Poisson_Solver_2D(xcc,zcc,-Jz0(:,:,n)',0,1);
        Az(:,:,n) = Az0';
        %
        display(n);
    end

else
    
    %%%   use Helmholtz theorem to calculate solonoidal and irrotational comps
    %%%   of electric field: E = -grad(phi) + curl(A)
    %
    theta = linspace(0,2*pi,20);
    dtheta = theta(2)-theta(1);
    theta = theta(1:end-1)+dtheta/2;
    Xmatrix = zeros(length(Zce),length(Xcc),length(theta));
    Ymatrix = zeros(length(Zce),length(Xcc),length(theta));
    Zmatrix = zeros(length(Zce),length(Xcc),length(theta));
    %
    XmatrixR = zeros(length(Zce),length(Xcc),length(theta));
    YmatrixR = zeros(length(Zce),length(Xcc),length(theta));
    %
    Ex_cez = zeros(size(divE_cez));
    Ez_cez = zeros(size(divE_cez));
    for j=1:length(Zce)
        for l=1:length(theta)
            Ymatrix(j,:,l) = Xcc*sin(theta(l));
            Xmatrix(j,:,l) = Xcc*cos(theta(l));
            %
            YmatrixR(j,:,l) = (Xce(end-1) + 0*Xcc)*sin(theta(l));
            XmatrixR(j,:,l) = (Xce(end-1) + 0*Xcc)*cos(theta(l));
        end
        Ex_cez(j,:,:) = (Ex(j+1,:,:)+Ex(j,:,:))/2.0;
        Ez_cez(j,:,:) = (Ez(j+1,:,:)+Ez(j,:,:))/2.0;
    end
    for i=1:length(Xcc)
        for l=1:length(theta)
            Zmatrix(:,i,l) = Zce;
        end
    end
    Rmatrix = sqrt(Xmatrix.^2+Ymatrix.^2);
    dV = Rmatrix*dX*dZ*dtheta;
    divE_3D = zeros(length(Zce),length(Xcc),length(tout),length(theta));
    curlE_3D = zeros(length(Zce),length(Xcc),length(tout),length(theta));
    Ex_3D   = zeros(length(Zce),length(Xcc),length(tout),length(theta));
    Ez_3D   = zeros(length(Zce),length(Xcc),length(tout),length(theta));
    for l=1:length(theta)
        divE_3D(:,:,:,l) = divE_cez;
        curlE_3D(:,:,:,l) = curlE_cez;
        Ex_3D(:,:,:,l)   = Ex_cez;    
        Ez_3D(:,:,:,l)   = Ez_cez; 
    end
    

    [~,ix0]=min(abs(Xcc-1));
    for i=ix0 %3:length(Xcc)-2
        for j=3:length(Zcc)-2
            for l=1:length(theta)
            rmag = sqrt( (Xcc(i)-Xmatrix).^2 + (Zcc(j)-Zmatrix).^2 ...
                 +       (0     -Ymatrix).^2 ); 
            rmagR = sqrt((Xcc(i)-XmatrixR).^2 + (Zcc(j)-Zmatrix).^2 ...
                  +      (0     -YmatrixR).^2); 
            rmag_pL = sqrt((Xcc(i)-Xmatrix).^2 + (Zcc(j)-Zce(end-1)).^2 ...
                    +      (0     -Ymatrix).^2); 
            rmag_mL = sqrt((Xcc(i)-Xmatrix).^2 + (Zcc(j)-Zce(2)).^2 ...
                    +      (0     -Ymatrix).^2);
            for k=itout0(end) %1:length(tout)
                int0 = squeeze(curlE_3D(:,:,k,:)).*dV./rmag;
                int1 = squeeze(divE_3D(:,:,k,:)).*dV./rmag;
                int2 = squeeze(Ex_3D(:,end-2,k,:)+Ex_3D(:,end-1,k,:))/2./squeeze(rmagR(:,i,:))*dtheta*Xce(end-1)*dZ;
                int3 = squeeze(Ez_3D(end-1,:,k,:)).*squeeze(Rmatrix(j,:,:)./rmag_pL(j,:,:))*dX*dtheta;
                int4 = squeeze(Ez_3D(2,:,k,:)).*squeeze(Rmatrix(j,:,:)./rmag_mL(j,:,:))*dX*dtheta;
                %
                thisSurfInt = sum(sum(int2(2:end-1,:)))/4/pi ...
                            + sum(sum(int3(3:end-2,:)-int4(3:end-2,:)))/4/pi;
                phi(j,i,k) = sum(sum(sum(int1(2:end-1,3:end-2,:))))/4/pi ...
                           - thisSurfInt;
                Cy(j,i,k)  = sum(sum(sum(int0(2:end-1,3:end-2,:))))/4/pi ...
                           - thisSurfInt;
            end
            end
            display(j);
        end
        display(i);
    end
    phi(1,:,:) = phi(end-3,:,:);
    phi(2,:,:) = phi(end-2,:,:);
    phi(end-1,:,:) = phi(3,:,:);
    phi(end,:,:) = phi(4,:,:);

end



%%%   compute elecrostatic field from phi
%
Ez_phi = zeros(size(Ez));
for j=2:nZ-1
   Ez_phi(j,:,:) = -(phi(j+1,:,:)-phi(j-1,:,:))/2/dZ;
end
Ez_phi(1,:,:) = Ez_phi(end-3,:,:);
Ez_phi(end,:,:) = Ez_phi(4,:,:);
Ez_em = Ez-Ez_phi;
%
Ex_phi = zeros(size(Ez));
for i=2:nX-1
   Ex_phi(:,i,:) = -(phi(:,i+1,:)-phi(:,i-1,:))/2/dX;
end
Ex_phi(1,:,:) = Ex_phi(end-3,:,:);
Ex_phi(end,:,:) = Ex_phi(4,:,:);
Ex_em = Ex-Ex_phi;

% figure(7); hold on; plot(Zcc/(2*pi)*3,phi0,'r'); grid on; box on;
% xlabel('z/\lambda'); title('electrostatic potential at r/a = 1');
% xlim([-0.5 0.5]); set(gca,'xtick',-0.5:0.25:0.5);
% %
% figure(8); hold on; plot(Zcc/(2*pi)*3,Ez(:,165,itout0(end))); grid on; box on;
% figure(8); hold on; plot(Zcc/(2*pi)*3,Ez_phi(:,165,itout0(end)),'r--'); 
% xlabel('z/\lambda'); title('z-electric field at r/a = 1');
% legend('total','electrostatic','location','best');
% xlim([-0.5 0.5]); set(gca,'xtick',-0.5:0.25:0.5);



%%%   plot potentials
%
f7=figure(7); set(f7,'position',[1200 700 1200 420]);
%
subplot(1,3,1);
h2=pcolor(Zcc,Xcc,phi(:,:,itout0(3))'); colorbar; box on
xlabel('z/a'); shading flat; ylabel('r/a'); 
title(['\phi at t = ',num2str(tout(itout0(3)),3),' t_0']);
axis('equal'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
colormap('jet');
%
subplot(1,3,2);
h2=pcolor(Zcc,Xcc,Ax(:,:,itout0(3))'); colorbar; box on
xlabel('z/a'); shading flat; ylabel('r/a'); 
title(['A_x at t = ',num2str(tout(itout0(3)),3),' t_0']);
axis('equal'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
colormap('jet');
%
subplot(1,3,3);
h2=pcolor(Zcc,Xcc,Az(:,:,itout0(3))'); colorbar; box on
xlabel('z/a'); shading flat; ylabel('r/a'); 
title(['A_z at t = ',num2str(tout(itout0(3)),3),' t_0']);
axis('equal'); axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
colormap('jet');




%%%   plot electric field 
%
close(figure(8));
f8=figure(8); set(f8,'position',[1200 700 1200 820]);
%f1=figure(1); set(f1,'position',[400 400 1800 200]);
set(gcf,'color','w');


subplot(2,3,1);
h1=pcolor(Zcc,Xcc,Ez_phi(:,:,itout0(end))'); colorbar; box on
xlabel('z/a'); shading flat; 
ylabel('r/a'); 
title(['-\partial\phi/\partial z at t = ',num2str(tout(itout0(end)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(2,3,2);
h2=pcolor(Zcc,Xcc,Ez_em(:,:,itout0(end))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['-\partial A_z/\partial t at t = ',num2str(tout(itout0(end)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(2,3,3);
h3=pcolor(Zcc,Xcc,Ez(:,:,itout0(end))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['E_z at t = ',num2str(tout(itout0(end)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
subplot(2,3,4);
h1=pcolor(Zcc,Xcc,Ex_phi(:,:,itout0(end))'); colorbar; box on
xlabel('z/a'); shading flat; 
ylabel('r/a'); 
title(['-\partial\phi/\partial x at t = ',num2str(tout(itout0(end)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(2,3,5);
h2=pcolor(Zcc,Xcc,Ex_em(:,:,itout0(end))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['-\partial A_x/\partial t at t = ',num2str(tout(itout0(end)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
subplot(2,3,6);
h3=pcolor(Zcc,Xcc,Ex(:,:,itout0(end))'); colorbar; box on
xlabel('z/a'); shading flat;
ylabel('r/a'); 
title(['E_x at t = ',num2str(tout(itout0(end)),3),' t_0']);
axis('equal'); 
axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
%
%
%
colormap('jet');    
%colormap('hot'); 


