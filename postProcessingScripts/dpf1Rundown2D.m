%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 1D dpf rundown module
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;




numProcs = 4;
%filePath1D = '../physicsMods/dpfRundown1D_2Temp/dataSave_1MAcyl/';
filePath1D = '../physicsMods/dpfRundown1D/resistiveMHD/1MAcyl_TVD2/'; numProcs = 4;
N0  = loadData1DVec(filePath1D,numProcs,'N');
M0  = loadData1DVec(filePath1D,numProcs,'M');
B0  = loadData1DVec(filePath1D,numProcs,'B');
P0  = loadData1DVec(filePath1D,numProcs,'P');
Ez0  = loadData1DVec(filePath1D,numProcs,'Ez');
Jz0 = loadData1DVec(filePath1D,numProcs,'J');
V0 = M0./N0;


numProcs = 50;
filePath = '../physicsMods/dpfRundown2D/resistiveMHD/testing1MAcyl/'; 
filePath = '../physicsMods/dpfRundown2D/resistiveMHD/testingCyl_C2/'; 
filePath = '../physicsMods/dpfRundown2D/resistiveMHD/1MACyl_TVD2/';

Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Zcc = loadData(filePath,numProcs,'Zcc');
Zce = loadData(filePath,numProcs,'Zce');
tout = loadData(filePath,numProcs,'tout');


N  = loadData(filePath,numProcs,'N');
Mx  = loadData(filePath,numProcs,'Mx');
By  = loadData(filePath,numProcs,'By');
P  = loadData(filePath,numProcs,'P');
Ez  = loadData(filePath,numProcs,'Ez');
Ex  = loadData(filePath,numProcs,'Ex');
Jz  = loadData(filePath,numProcs,'Jz');
Jx  = loadData(filePath,numProcs,'Jx');
Vx  = loadData(filePath,numProcs,'Vx');
Vz  = loadData(filePath,numProcs,'Vz');
Jxcc  = loadData(filePath,numProcs,'Jxcc');
Jzcc  = loadData(filePath,numProcs,'Jzcc');
eta  = loadData(filePath,numProcs,'eta');
hy_cc  = loadData(filePath,numProcs,'hy_cc');
hy_ce  = loadData(filePath,numProcs,'hy_ce');
%V = Mx./N;
T = P./N/2.0;



%%%   plot contours
%
close(figure(1));
f1=figure(1); set(f1,'position',[1100 300 1500 500]); % for ka=3
%f1=figure(1); set(f1,'position',[1500 30 1080 770]); % for ka=10
%f1=figure(1); set(f1,'position',[400 400 1800 200]);
set(gcf,'color','w');

thisFigNum = 0;
for thist=1:length(tout)

    subplot(1,3,1);
    h1=pcolor(Zce(2:end-1),Xce(2:end-1),N(3:end-1,3:end-1,thist)'); colorbar; box on
    xlabel('z/a'); shading flat;
    ylabel('r/a'); 
    title(['density at t = ',num2str(tout(thist),3),' t_0']);
    axis('equal'); 
    axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
    %
    %
    %
    subplot(1,3,2);
    h2=pcolor(Zce(2:end-1),Xce(2:end-1),By(3:end-1,3:end-1,thist)'); colorbar; box on
    xlabel('z/a'); shading flat;
    ylabel('r/a'); 
    title(['B_y at t = ',num2str(tout(thist),3),' t_0']);
    axis('equal'); 
    axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
    colormap('jet');
    %
    %
    %
    subplot(1,3,3);
    h3=pcolor(Zce(2:end-1),Xce,Ez(3:end-1,:,thist)'); colorbar; box on
    xlabel('z/a'); shading flat;
    ylabel('r/a'); 
    title(['E_z at t = ',num2str(tout(thist),3),' t_0']);
    axis('equal'); 
    axis([Zce(2) Zce(end-1) Xce(3) Xce(end-2)]);
    colormap('jet');    
    
%     %%%   put time stamp on figure
%     %
%     a1=annotation(f1,'textbox',...
%     [0.03 0.5 0.1 0.05],...
%     'String',['t = ',num2str(tout(thist)),' t_A'],...
%     'FitBoxToText','off', ...
%     'backgroundcolor','y');
    
    
    M(thisFigNum+1) = getframe(f1); 
    thisFigNum = thisFigNum+1;
    if(tout(thist)~=tout(end))
      %  delete(a1);
        delete(h1);
        delete(h2);
        delete(h3);
    end

end

% f1pos = get(f1,'position');
% f11=figure(11);
% set(f11,'Position',f1pos);
% movie(f11,M,1,3);
% %v=VideoWriter([dataPath,'movieFigs/2DrBthMovie.avi']);
% v=VideoWriter('./2Dpinch.avi');
% v.FrameRate = 4; %v.Quality=100; %v.CompressionRatio = 2;
% %v.LosslessCompression = true;
% open(v); writeVideo(v,M);
% close(v);








%%%   compute z-symmetric error for density
%
errorN = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
errorVx = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
errorVz = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
errorJxcc = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
%errorEx = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
errorP = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
errorBy = zeros((length(Zcc)-4)/2,length(Xcc),length(tout));
for j=1:length(errorN(:,1,1))
    errorN(j,:,:) = N(j+2,:,:)  - N(end-1-j,:,:);
    errorVx(j,:,:) = Vx(j+2,:,:)  - Vx(end-1-j,:,:);
    errorVz(j,:,:) = Vz(j+2,:,:)  + Vz(end-1-j,:,:);
    errorJxcc(j,:,:) = Jxcc(j+2,:,:)  + Jxcc(end-1-j,:,:);
    errorBy(j,:,:) = By(j+2,:,:) - By(end-1-j,:,:);
   % errorEx(j,:,:) = S(j+2,:,:)  - S(end-1-j,:,:);
    errorP(j,:,:) = P(j+2,:,:)  - P(end-1-j,:,:);
end

figure(7); 
pcolor(Zcc(3:length(errorN(:,1,1))+2),Xcc,log10(abs(errorN(:,:,end)))'); 
shading flat; colorbar; title('error in N symmetry');

figure(8); 
pcolor(Zcc(3:length(errorN(:,1,1))+2),Xcc,log10(abs(errorBy(:,:,end)))'); 
shading flat; colorbar; title('error in P symmetry');

colormap('jet');

maxErrorN0 = max(max(errorN(:,:,1)))
maxErrorP0 = max(max(errorP(:,:,1)))
