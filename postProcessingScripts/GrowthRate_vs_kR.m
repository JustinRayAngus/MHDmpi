%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   plot seeded m=0 mode growth rate for different kR
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


filePath = '../physicsMods/pinch2D/';
newDeck=0;


%%%   kR = 1
%
%filePath = '../../fromQuartz/pinch2D/kR1/data_80nz_400nx/'; numProcs = 20;
filePath = '../../fromQuartz/pinch2D/kR1/data_160nz_800nx/'; numProcs = 80;
%filePath = '../../fromQuartz/pinch2D/kR1/data_320nz_800nx/'; numProcs = 80;

%%%   kR = 5
%
filePath = '../../fromQuartz/pinch2D/kR5/data_160nz_800nx/'; numProcs = 80;


%%%   kR = 10
%
filePath = '../../fromQuartz/pinch2D/kR10/data_160nz_800nx/'; numProcs = 80;


%%%   kR = 20
%
%filePath = '../../fromQuartz/pinch2D/kR20/data_80nz_400nx/'; numProcs = 20;


%%%   kR = 30
%
%filePath = '../../fromQuartz/pinch2D/kR30/data_80nz_400nx/'; numProcs = 20;


%%%   kR = 40
%
%filePath = '../../fromQuartz/pinch2D/kR40/data_80nz_400nx/'; numProcs = 20;


%%%   Hall MHD (hall plus grad Pe)
%
%filePath = '../../fromQuartz/pinch2D/Hall_v2/kR10/test5_Li=5.0e-2/'; numProcs = 20;
%filePath = '../../fromQuartz/pinch2D/Hall_v2/kR10/test5_Li=5.0e-3/'; numProcs = 20;


%%%   drift-ideal MHD
%
filePath = '../../fromQuartz/pinch2D/entropy_v0/Li0.0/ka20/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v0/Li1.5e-2/ka20/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v0/Li1.5e-2/ka3.0_stable_nuT/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v0/Li5.0e-2/ka20_smallerDt/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v0/Li5.0e-2/ka3.0/'; numProcs = 20; newDeck = 2;


filePath = '../../fromQuartz/pinch2D/entropy_v0/Li0.12/ka10/'; numProcs = 20; newDeck = 2;
%filePath = '../../fromQuartz/pinch2D/driftIdeal_v2/ka3.0/'; numProcs = 20; newDeck = 2;
filePath = '../../fromQuartz/pinch2D/entropy_v1/Li0.12/noGyroVisc/ka10.0_taui1.0e-2/'; numProcs = 20; newDeck = 2;


t0 = 1.2046e-8; % see normalizationParameters_zPinch.m
tA = 2.5e-8;    % from PoP 17, 072107 (2010)
%tA = 5.4036e-08;  % r0/VTi
%tA = 3.8209e-08;


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
Jx  = loadData(filePath,numProcs,'Jx');
Ez  = loadData(filePath,numProcs,'Ez');
Ex  = loadData(filePath,numProcs,'Ex');
Cs  = loadData(filePath,numProcs,'Cs');
gamma0 = loadData(filePath,numProcs,'gamma0');
dX = Xcc(2)-Xcc(1);
dZ = Zcc(2)-Zcc(1);


%%%   calculate global values for conservation
%
dV = rcc*2*pi*dX*dZ;
Mass = zeros(size(tout));
Entropy = zeros(size(tout));
zMomentum = zeros(size(tout));
ByFlux = zeros(size(tout));
Energy = zeros(size(tout));
Edens = 0.5*(Mx.*Vx+Mz.*Vz) + 1.5*P + By.^2/2.0; 

for n=1:length(tout)
    NdV  = N(:,:,n).*dV;
    SdV  = S(:,:,n).*dV;
    MzdV = Mz(:,:,n).*dV;
    EdensdV = Edens.*dV;
    Mass(n) = sum(sum(NdV(3:end-2,3:end-2)));
    Entropy(n) = sum(sum(SdV(3:end-2,3:end-2)));
    ByFlux(n) = sum(sum(By(3:end-2,3:end-2,n)))*dX*dZ;
    zMomentum(n) = sum(sum(MzdV(3:end-2,3:end-2)));
    Energy(n) = sum(sum(EdensdV(3:end-2,3:end-2)));
end


close(figure(4));
figure(4); plot(tout,Mass,'displayName','mass');
hold on; plot(tout,Entropy,'displayName','entropy');
hold on; plot(tout,zMomentum,'displayName','z-momentum');
hold on; plot(tout, ByFlux,'displayName','B_y-flux');
xlabel('t/t_0'); title('global conservations');
lg4 = legend('show');
set(lg4,'location','best');


%%%   calculate perturbed mass density per unit length 
%
dA = rcc*2*pi*dX;
deltaN = zeros(size(N));
Dz = zeros(length(Zcc)-4,length(tout)); % mass per unit length z


for n=1:length(tout)
    deltaN(:,:,n) = N(:,:,n);
 %   deltaN(:,:,n) = Te(:,:,n);
    %  deltaN(:,:,n) = S(:,:,n);
    deltaNdA = squeeze(deltaN(3:end-2,3:end-2,n).*dA(3:end-2,3:end-2));
    Dz(:,n) = sum(deltaNdA,2);
end
%Dz = squeeze(N(3:end-2,60,:));


%%%   calculate Fourier transform of all modes and amplitudes
%%%   assumed period on Zce grid
%
Nz = length(Zce)-2;
L = Zce(end-1)-Zce(2);

Zmodes = 0:1:(Nz-1)/2; % assuming N is odd
k = 2*pi*Zmodes/L;
DzFT = zeros(length(k),length(tout)); % Fourier transform

for n=1:length(tout)
   for j=1:length(Zmodes)
      Dzexp = Dz(:,n).*exp(-1i*k(j)*Zcc(3:end-2));
      DzFT(j,n) = sum(Dzexp)*dZ;
   end
end
DzFT_amp = abs(DzFT)/Mass(1);
%DzFT_amp = abs(real(DzFT))/Mass(1);
DzFT_amp_odd = abs(imag(DzFT))/Mass(1);

%tout = sqrt(2)*tout;
f1=figure(2); set(f1,'position',[1450 600 1000 450]);
%
subplot(1,2,1);
hold on; plot(tout,DzFT_amp(2,:),'displayName',['kR=',num2str(k(2))]);
set(gca,'yscale','log'); box on;
%hold on; plot(tout,DzFT_amp(3,:),'displayName',['kR=',num2str(k(3))]);
%hold on; plot(tout,DzFT_amp(4,:),'displayName',['kR=',num2str(k(4))]);
xlabel('t/t_0'); ylabel('mode amplitude'); 
title('Fourier Amplitudes'); grid on, grid off; grid on;
lg9=legend('show'); set(lg9,'location','best');
axis([0 tout(end) 1.0e-6 1]); axis('square');
set(gca,'ytick',10.^(-10:2:0));



%%% for looking at modes of N
%hold on; plot(tout,2.8e-5*10.^(tout/2.15),'b--');
%hold on; plot(tout,2.0e-10*10.^(tout*2/2.15),'r--');
%hold on; plot(tout,7.5e-19*10.^(tout*3/2.15),'linestyle','--'); ? should
%be 3 times as fast, not 4 times as fast


%%% for looking at modes of S
%
%hold on; plot(tout,4.2e-5*10.^(tout/2.15),'b--');
%hold on; plot(tout,8.5e-10*10.^(tout*2/2.15),'r--');
%hold on; plot(tout,2.0e-14*10.^(tout*3/2.15),'linestyle','--');
%hold on; plot(tout,0.9e-5*10.^(tout*0.06)); % Hall_v2, kR=30

%%%   calculate slope of Fourier modes
%
mDzFT_amp = zeros(length(k),length(tout)-1);
for i=1:length(k)
    for j=1:length(tout)-1
        mDzFT_amp(i,j) = (log(DzFT_amp(i,j+1))-log(DzFT_amp(i,j)))/(tout(j+1)-tout(j));
    end
end



subplot(1,2,2);
hold on; plot(tout(1:end-1),mDzFT_amp(2,:),'displayName',['kR=',num2str(k(2))]); 
box on;
xlabel('t/t_0'); ylabel('growth rate \times t_0');
title('derivative of log Fourier Amps');
%axis([0 10 0 2]); 
axis('square'); grid on;



NperZ = Mass(1)/dZ/length(Zcc(3:end-2))/pi; % ave density per unit length
Mi = 1.008*1.6605402e-27; % [kg]
n0 = 1.0e24; % [1/m^3]
R0 = 5.0e-3; % [m]

rhobar = Mi*NperZ*n0; % [kg/m^3]


%%%   t0 = a/Cs, with Cs=sqrt(2*T0/Mi)   
%
ka = [0.0 0.3 1 3 5 10 15 20];
gamma_Li0      = [0 0.123 0.33 0.50 0.54 0.57 0.57 0.57]; % Li/a = 0.0;
gamma_Li1p5em2 = [0 0.123 0.33 0.51 0.56 0.62 0.66 0.71]; % Li/a = 1.5e-2
gamma_Li5p0em2 = [0 0.127      0.57                1.00]; % Li/a = 5.0e-2
gamma_Li0p12   = [0 0.136 0.40 0.68 0    1.00 0    0   ]; % Li/a = 0.12
kasub = ka([1,3:4,6]);


%%%   values form Kurt (LSP) for FUZE (Li/a=0.12)
%%%   t0 = a/Va, with Va calculated using peak B (t0=2.6ns)
kaKurt = [0 1 1.5 2 2.5 5 7.5 10.0 12.5 15];
gamma_Kurt = [0 0.49 0.63 0.68 0.70 0.78 0.71 0.59 0.51 0.47]/sqrt(2);
%
kaVasi = [2.5 5 7.5 10];
gamma_Vasi = [1.0 1.05 0.9 0.65]/sqrt(2);

close(figure(3)); f3=figure(3); 
plot(ka,gamma_Li0,'Marker','*');grid on;
hold on; plot(ka,0.6+0*ka); % max gamma from local analysis
xlabel('ka'); ylabel('\gammat_0'); title('Bennett Equilibrium Growth Rates');
axis('square'); axis([0 20 0 0.7]);




