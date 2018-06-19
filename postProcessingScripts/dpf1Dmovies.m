%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   This script makes movies of 1D dpf rundown sims
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
addpath('~angus1/Programs/MATLAB_LSP_TOOLS/');


numProcs = 4;
filePath = '../physicsMods/dpfRundown1D/';


%%%   set default font and lines
%
set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'defaultaxeslinewidth',1.5);
set(0,'defaultlinelinewidth',3);
set(0,'defaultaxesfontweight','bold');

deleteFlag = 0;
thisFigNum = 0;
for it = 1:201

    f1=figure(1); 
    set(f1,'position',[1030 425 1200 500]);
    set(gcf,'color','w');
    
    for i=1:numProcs
        fileName = ['output',num2str(i-1),'.h5'];
        thisFile = [filePath,fileName];
        procID  = hdf5read(thisFile,'procID');
        fileinfo = hdf5info(thisFile);
        Xcc = hdf5read(thisFile,'Xcc');
        Xce = hdf5read(thisFile,'Xce');
        N  = hdf5read(thisFile,'N');
        M  = hdf5read(thisFile,'M');
        S  = hdf5read(thisFile,'S');
        B  = hdf5read(thisFile,'B');
        P  = hdf5read(thisFile,'P');
        V  = hdf5read(thisFile,'V');
        J  = hdf5read(thisFile,'J');
        J0  = hdf5read(thisFile,'J0'); % curl(B)
        Ez  = hdf5read(thisFile,'Ez');
        Cs  = hdf5read(thisFile,'Cs');
        eta  = hdf5read(thisFile,'eta');
        gamma0 = hdf5read(thisFile,'gamma0');
        FluxRatio  = hdf5read(thisFile,'FluxRatio');
        FluxLim    = hdf5read(thisFile,'FluxLim');
        FluxL    = hdf5read(thisFile,'FluxL');
        FluxR    = hdf5read(thisFile,'FluxR');
        FluxN  = hdf5read(thisFile,'FluxN');
        FluxM  = hdf5read(thisFile,'FluxM');
        FluxS  = hdf5read(thisFile,'FluxS');
        FluxB  = hdf5read(thisFile,'FluxB');
        FluxEz  = hdf5read(thisFile,'FluxEz');
        tout= hdf5read(thisFile,'tout');
        dX = Xcc(2)-Xcc(1);
        T = P./N;
        Ptot = P+B.^2/2+M.^2./N;
        
        subplot(1,3,1);
        hold on; plot(Xcc,N(:,it),'black'); box on; grid on;
        set(gca,'xtick',0:0.1:0.5);
        %set(gca,'ytick',0:0.3:1.2);
        xlabel('x'); ylabel('N');
        title('mass density'); axis('square');
        axis([0 0.5 0 4]);
        %
        subplot(1,3,2);
        hold on; plot(Xcc,Ptot(:,it),'black'); box on; grid on;
        hold on; plot(Xcc,B(:,it).^2/2,'b');
        hold on; plot(Xcc,P(:,it),'r');
        hold on; plot(Xcc,M(:,it).^2./N(:,it),'color',[0.47 0.67 0.19]);
        set(gca,'xtick',0:0.1:0.5);
        %set(gca,'ytick',0:0.3:1.2);
        xlabel('x'); ylabel('P');
        title('pressure'); axis('square');
        if(i==1)
            legend('total','magnetic','thermal','mean','location','NW');
        end
        axis([0 0.5 0 10]);
        %
        %
        %
        subplot(1,3,3);
        hold on; plot(Xce,J(:,it),'black'); box on; grid on;
        xlabel('x'); ylabel('J'); axis('square');
        title('current density');
        set(gca,'xtick',0:0.1:0.5);
        %set(gca,'ytick',1:0.5:3);
        axis([0 0.5 -25 0]);

    end

    Mov(thisFigNum+1) = getframe(f1); 
 %   filename = ['fig',num2str(thisFigNum,'%03d')];
 %   saveas(f1,[savepath,filename],'png');
    thisFigNum = thisFigNum+1;
    close(f1);
    
end


%close(figure(11));
f2=figure(2);
set(f2,'Position',[1030 425 1200 500]);
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


