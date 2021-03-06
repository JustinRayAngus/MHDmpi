%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   This script makes movies of 1D dpf rundown sims
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
addpath('~angus1/Programs/MATLAB_LSP_TOOLS/');


%%%   set default font and lines
%
set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'defaultaxeslinewidth',1.5);
set(0,'defaultlinelinewidth',3);
set(0,'defaultaxesfontweight','bold');


%%%   define number of procs and load date
%
numProcs = 4;
%filePath = '../physicsMods/dpfRundown1D_Econs/';
filePath = '../physicsMods/dpfRundown1D_2Temp/dataSave_1MAcar/'; TwoTempVersion=1;
filePath = '../physicsMods/dpfRundown1D_2Temp/dataSave_1MAcyl/'; TwoTempVersion=1;
%filePath = '../physicsMods/dpfRundown1D_2Temp/';


Xcc = loadData1DVec(filePath,numProcs,'Xcc');
Xce = loadData1DVec(filePath,numProcs,'Xce');
tout = loadData1DVec(filePath,numProcs,'tout');

% fileName = ['output',num2str(i-1),'.h5'];
% thisFile = [filePath,fileName];
% procID  = hdf5read(thisFile,'procID');
% fileinfo = hdf5info(thisFile);

Btot = 0*tout;
Poynt = 0*tout;
Bflux = 0*tout;
Mflux = 0*tout;
Msource = 0*tout;
Ssource = 0*tout;
Psource = 0*tout;
Poynt_exit = 0*tout;
EEztot = 0*tout;

N  = loadData1DVec(filePath,numProcs,'N');
M  = loadData1DVec(filePath,numProcs,'M');
B  = loadData1DVec(filePath,numProcs,'B');
P  = loadData1DVec(filePath,numProcs,'P');
if(TwoTempVersion)
    Pe  = loadData1DVec(filePath,numProcs,'Pe');
    Pi  = loadData1DVec(filePath,numProcs,'Pi');
    Te  = loadData1DVec(filePath,numProcs,'Te');
    Ti  = loadData1DVec(filePath,numProcs,'Ti');
    %
    hy_cc = loadData1DVec(filePath,numProcs,'hy_cc');
    hy_ce = loadData1DVec(filePath,numProcs,'hy_ce');
else
    T  = loadData1DVec(filePath,numProcs,'T');
    Te = T/2; % units are 1/2 eV
    Ti = T/2; % units are 1/2 eV
    Pe = P/2;
    Pi = P/2;
    %
    hy_cc = 1.0 + 0.0*Xcc;
    hy_ce = 1.0 + 0.0*Xce;
end
Te = 2.0*Te; % [eV]
Ti = 2.0*Ti; % [eV]
T = (Te+Ti)/2.0;

V  = loadData1DVec(filePath,numProcs,'V');
J  = loadData1DVec(filePath,numProcs,'J');
Jcc = loadData1DVec(filePath,numProcs,'Jcc');
J0  = loadData1DVec(filePath,numProcs,'J0');
Ez  = loadData1DVec(filePath,numProcs,'Ez');
eta = loadData1DVec(filePath,numProcs,'eta');
gamma0 = loadData1DVec(filePath,numProcs,'gamma0');
delta0 = loadData1DVec(filePath,numProcs,'delta0');
FluxM = loadData1DVec(filePath,numProcs,'FluxM');
FluxN = loadData1DVec(filePath,numProcs,'FluxN');

%T = P/2.0./N; % average temperature

%%%   calculate divU
%
dX = Xcc(2)-Xcc(1);
dV = dX*hy_cc;
divV = zeros(size(V));
dPdx = zeros(size(V));
dhydx = zeros(size(Xcc));
for j=2:length(Xcc)-1
    divV(j,:) = (V(j+1,:)-V(j-1,:))/2/dX;
    dPdx(j,:) = (P(j+1,:)-P(j-1,:))/2/dX;
    dhydx(j) = (hy_ce(j)-hy_ce(j-1))/dX;
end
dhydx(1) = dhydx(2);
dhydx(end) = dhydx(end-1);


%%%   calcualte Ez at cell-center
%
Ezcc = zeros(size(N));
for j=3:length(Xcc)-2
    Ezcc(j,:) = (Ez(j-1,:) + Ez(j,:))/2;
end


%figure(77); hold on; plot(Xcc,divV(:,100)); box on;

Pmag = B.^2/2;
Pmageff = Pmag + cumtrapz(Xcc,B.^2./Xcc); % for cyl only
Pram = M.*V;
Ptot = P + B.^2/2 + M.*V/2.0;
PdV = P.*divV;
VdP = V.*dPdx;
etaJ2 = eta.*Jcc.^2;
etaJ2mod = etaJ2./T;
Edenstot = 3/2*P + N.*V.^2/2 + B.^2/2;

%%%   calculate conservation stuff
%
JdotE  = sum(Jcc(3:end-2,:).*Ezcc(3:end-2,:).*dV(3:end-2));
Etot   = sum(Edenstot(3:end-2,:).*dV(3:end-2));
Etherm = sum(1.5*P(3:end-2,:).*dV(3:end-2));
Emean  = sum(0.5*M(3:end-2,:).*V(3:end-2,:).*dV(3:end-2));
Efield = sum(B(3:end-2,:).^2/2.*dV(3:end-2));
EEztot = delta0*sum(Ez(2:end-2,:).^2/2.0.*dV(3:end-2));
%

EntropyDensity = 2.0*N.*log(T.^1.5./N);
Entropy  = sum(EntropyDensity(3:end-2,:).*dV(3:end-2));
Momentum = sum(M(3:end-2,:).*dV(3:end-2));
Mass     = sum(N(3:end-2,:).*dV(3:end-2));

for n=1:length(tout)

    Btot(n)   = sum(B(3:end-2,n))*dX;
    Ssource(n) = sum(etaJ2mod(3:end-2,n).*dV(3:end-2));
    Msource(n) = sum(-Jcc(3:end-2,n).*B(3:end-2,n).*dV(3:end-2) ...
                     +P(3:end-2,n).*dhydx(3:end-2)*dX);
    Psource(n) = sum(etaJ2(3:end-2,n)-0*PdV(3:end-2,n).*dV(3:end-2));
   % Psource(n) = sum(etaJ2(3:end-2,n)+VdP(3:end-2,n))*dX;
    Etot(n)    = Etot(n) + EEztot(n);
    
    Poynt(n) = Ez(end-1,n)*(B(end-2,n)+B(end-1,n))/2.0 ...
             - Ez(2,n)*(B(3,n)+B(2,n))/2.0 ...
             - 0*5/2*(P(2,n)+P(3,n))/2*(V(2,n)+V(3,n))/2 ...
             - 0*0.5*N(3,n)*V(3,n)^3;
    Bflux(n) = Ez(end-1,n) - Ez(2,n);
    Mflux(n) =  (P(2,n)+P(3,n))/2 + B(2,n)^2/2;
    Msource(n) = Msource(n) + FluxM(2,n) - FluxM(end-1,n);
  %  Psource(n) = Psource(n) - 3/2*P(2,n)*(V(2,n)+V(3,n))/2;

    Msource(n) = Msource(n) - (P(end-2,n)+P(end-1,n))/2;
    Mflux(n) = Mflux(n) - (P(end-2,n)+P(end-1,n))/2;
   
end
EthermGain = cumtrapz(tout,Psource);
Cs = sqrt(gamma0*P./N);
Mach = abs(V)./Cs;


% rho2/rho1 = v1/v2 = (gamma0+1)/(gamma0-1+2/M1^2)
M1 = 1.1;
deltaN = (gamma0+1)/(gamma0-1+2/M1^2);
deltaP = (2*gamma0*M1.^2-(gamma0-1))/(gamma0+1);

Ez0 = eta.*Jcc-V.*B;

%%%%%%%%%%%%
%%%

deleteFlag = 0;
thisFigNum = 0;
itmax = 111; %length(tout);
for it = 1:1:itmax
    
f1=figure(11); 
set(f1,'position',[450 80 1340 440]);
set(gcf,'color','w');
    
subplot(1,3,1);
plot(Xcc,N(:,it),'black'); box on; grid on;
set(gca,'xtick',0:0.2:1);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x/R'); ylabel('N');
title('density'); axis('square');
%axis([0 1 0 12]);
axis([0 1 0 80]);


subplot(1,3,2);
%h4=plot(Xcc,B(:,it).^2/2+P(:,it),'black'); box on; grid on;
h4=plot(Xcc,B(:,it).^2/2,'black');
hold on; plot(Xcc,P(:,it),'r');
set(gca,'xtick',0:0.2:1);
%set(gca,'ytick',0:0.3:1.2);
xlabel('x/R'); ylabel('P');
title('pressure'); axis('square');
legend('magnetic','thermal','location','northeast');
grid on;
%axis([0 1 0 2000]);
axis([0 1 0 8e4]);


subplot(1,3,3);
plot(tout(1:it),Efield(1:it),'displayName','magnetic','color','black');
hold on; plot(tout(1:it),Etherm(1:it)-Etherm(1),'displayName','thermal','color','r'); 
hold on; plot(tout(1:it),Emean(1:it),'displayName','mean', ...
                   'color','b','linestyle','--');
lg3=legend('show','location','northwest');
xlabel('time'); ylabel('energy');
title('energy partition'); grid on;
axis('square');
set(gca,'xtick',0:0.02:0.1);
%axis([0 tout(itmax) 0 300]);
axis([0 tout(itmax) 0 700]);

%
%
%


%%%   put time stamp on figure
%
a1=annotation(f1,'textbox',...
[0.914 0.76 0.06 0.0743],...
'String',['t=',num2str(tout(it),3)],...
'FitBoxToText','off', ...
'backgroundcolor','y');


Mov(thisFigNum+1) = getframe(f1); 
%   filename = ['fig',num2str(thisFigNum,'%03d')];
%   saveas(f1,[savepath,filename],'png');
thisFigNum = thisFigNum+1;
if(it~=itmax)
    close(f1);
end
    
end


%close(figure(11));
f2=figure(2);
f1pos = get(f1,'position');
set(f2,'Position',f1pos);
movie(f2,Mov,1,3);
%v=VideoWriter([dataPath,'movieFigs/2DrBthMovie.avi']);
v=VideoWriter('./dpf1Dmovie.avi');
v.FrameRate = 4; %v.Quality=100; %v.CompressionRatio = 2;
%v.LosslessCompression = true;
open(v); writeVideo(v,Mov);
close(v);
%movie2avi(M,'movieTesting2.avi','fps',5);

%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


