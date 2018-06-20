%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   This script makes movies of 1D dpf rundown sims in cyl coords
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
addpath('~angus1/Programs/MATLAB_LSP_TOOLS/');


% %%%   set default font and lines
% %
% set(0,'defaultaxesfontsize',18);
% set(0,'defaulttextfontsize',18);
% set(0,'defaultaxeslinewidth',1.5);
% set(0,'defaultlinelinewidth',3);
% set(0,'defaultaxesfontweight','bold');


%%%   define number of procs and load date
%
numProcs = 4;
filePath = '../physicsMods/dpfRundown1Dcyl/';

Xcc = loadData1DVec(filePath,numProcs,'Xcc');
Xce = loadData1DVec(filePath,numProcs,'Xce');
tout = loadData1DVec(filePath,numProcs,'tout');
%
N  = loadData1DVec(filePath,numProcs,'N');
M  = loadData1DVec(filePath,numProcs,'M');
S  = loadData1DVec(filePath,numProcs,'S');
B  = loadData1DVec(filePath,numProcs,'B');
P  = loadData1DVec(filePath,numProcs,'P');
V  = loadData1DVec(filePath,numProcs,'V');
J  = loadData1DVec(filePath,numProcs,'J');
Jcc = loadData1DVec(filePath,numProcs,'Jcc');
J0  = loadData1DVec(filePath,numProcs,'J0');
Ez  = loadData1DVec(filePath,numProcs,'Ez');
eta = loadData1DVec(filePath,numProcs,'eta');
gamma0 = loadData1DVec(filePath,numProcs,'gamma0');
delta0 = loadData1DVec(filePath,numProcs,'delta0');
dX = Xcc(2)-Xcc(1);
%
%%%   calcualte Ez at cell-center
%
Ezcc = zeros(size(N));
for j=3:length(Xcc)-2
    Ezcc(j,:) = (Ez(j-1,:) + Ez(j,:))/2;
end
%
dV     = 2*pi*Xcc(3:end-2)*dX;
JdotE  = sum(Jcc(3:end-2,:).*Ezcc(3:end-2,:).*dV);
Etherm = sum(1.5*P(3:end-2,:).*dV);
Emean  = sum(0.5*M(3:end-2,:).*V(3:end-2,:).*dV);
Jheating = cumtrapz(tout,JdotE);
%
T = P./N;
Ptot = P+B.^2/2+M.^2./N;

Cs = sqrt(gamma0*P./N);
Mach = abs(V./Cs);

%%% calculate 


deleteFlag = 0;
thisFigNum = 0;
for it = 1:length(tout)
    
f1=figure(1); 
set(f1,'position',[450 80 1340 730]);
set(gcf,'color','w');
    
subplot(2,3,1);
plot(Xcc,N(:,it),'black'); box on; grid on;
set(gca,'xtick',0:0.2:1);
%set(gca,'ytick',0:0.3:1.2);
xlabel('r'); ylabel('N');
title('mass density'); axis('square');
axis([0 1 0 20]); axis('square');
%
subplot(2,3,2);
h4=plot(Xcc,B(:,it).^2/2+P(:,it),'black'); box on; grid on;
hold on; plot(Xcc,B(:,it).^2/2,'b');
hold on; plot(Xcc,P(:,it),'r');
set(gca,'xtick',0:0.2:1);
%set(gca,'ytick',0:0.3:1.2);
xlabel('r'); ylabel('P');
title('pressure'); axis('square');
legend('total','magnetic','thermal','location','best');
axis([0 1 0 80]);  axis('square');
%
%
%
subplot(2,3,3);
h3=plot(Xcc,Mach(:,it),'black'); box on; grid on;
xlabel('r'); ylabel('V/C_s'); axis('square');
title('local Mach Number');
set(gca,'xtick',0:0.2:1);
%set(gca,'ytick',1:0.5:3);
axis([0 1 0 1.5]);  axis('square');
%
%
%
subplot(2,3,4);
h4=plot(tout(1:it),Jheating(1:it),'displayName','\int \int J\cdot E dx dt'); box on; grid on;
hold on; plot(tout(1:it),Emean(1:it)+Etherm(1:it)-Etherm(1),'linestyle','--','displayName','\int (3/2P + \rhoV^2/2) dx');
xlabel('time'); ylabel('energy'); axis('square');
title('kinetic energy conservation');
%set(gca,'xtick',0:0.1:0.2);
%set(gca,'ytick',1:0.5:3);
axis([0 tout(end) 0 60]);  axis('square');
lg4 = legend('show','location','best');


%%%   put time stamp on figure
%
a1=annotation(f1,'textbox',...
[0.035 0.5 0.06 0.047],...
'String',['t=',num2str(tout(it),3)],...
'FitBoxToText','off', ...
'backgroundcolor','y');


Mov(thisFigNum+1) = getframe(f1); 
%   filename = ['fig',num2str(thisFigNum,'%03d')];
%   saveas(f1,[savepath,filename],'png');
thisFigNum = thisFigNum+1;
if(it~=length(tout))
close(f1);
end
    
end


%close(figure(11));
f2=figure(2);
f1pos = get(f1,'position');
set(f2,'Position',f1pos);
movie(f2,Mov,2,3);
%v=VideoWriter([dataPath,'movieFigs/2DrBthMovie.avi']);
v=VideoWriter('./dpf1Dmovie.avi');
v.FrameRate = 2; %v.Quality=100; %v.CompressionRatio = 2;
%v.LosslessCompression = true;
open(v); writeVideo(v,Mov);
close(v);
%movie2avi(M,'movieTesting2.avi','fps',5);

%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


