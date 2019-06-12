%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   look at output variables from 2D liftOff sims in rth
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%        create physical grid
%%%
%%%
Xcc = loadData(filePath,numProcs,'Xcc');
Xce = loadData(filePath,numProcs,'Xce');
Ycc = loadData(filePath,numProcs,'Zcc');
Yce = loadData(filePath,numProcs,'Zce');

r  = Xce(2:end-1);
th = Yce(2:end-1);

nr = length(r);
nth = length(th);

x = zeros(nr,nth); y = zeros(nr,nth);
xi0 = zeros(1,nth); yi0 = zeros(1,nth);
xi1 = zeros(1,nth); yi1 = zeros(1,nth);
ri0 = 0.8; ri1 = 1.0;
for j=1:nth
    for i=1:nr
        x(i,j) = r(i)*cos(th(j));
        y(i,j) = r(i)*sin(th(j));
        %    
    end
    xi0(j) = ri0*cos(th(j));
    yi0(j) = ri0*sin(th(j));    
    xi1(j) = ri1*cos(th(j));
    yi1(j) = ri1*sin(th(j));   
end
xbox_ins = [ri0 xi1 flip(xi0)];
ybox_ins = [0   yi1 flip(yi0)];


%%%   convert Bx ,By, Ez, to cell-center
%
nYcc = length(Ycc);
nXcc = length(Xcc);
Bx_cc = zeros(nr,nth,length(tns));
By_cc = zeros(nr,nth,length(tns));
for i=1:nr-1
    Bx_cc(i,1:nth-1,:) = (Bx(3:nYcc-2,2+i,:)+Bx(3:nYcc-2,3+i,:))/2;
end
Bx_cc(end,:,:) = Bx_cc(end-1,:,:);
Bx_cc(:,end,:) = Bx_cc(:,end-1,:);
for j=1:nth-1
    By_cc(1:nr-1,j,:) = (By(2+j,3:nXcc-2,:)+By(3+j,3:nXcc-2,:))/2;
end
By_cc(end,:,:) = By_cc(end-1,:,:);
By_cc(:,end,:) = By_cc(:,end-1,:);

close(figure(1));
f1=figure(1);
%set(f1,'Position',[1600 800 800 400]);
set(f1,'Position',[1500 800 1000 350]);
set(gcf,'color','w');

thisFigNum = 0;
for it = 1:length(tns)
    
    thist = tns(it);%

    subplot(1,2,1);
    minVal = min(min(Bx_cc(:,:,it)));
    maxVal = max(max(Bx_cc(:,:,it)));
    p1=pcolor(x,y,Bx_cc(:,:,it)); shading flat; cb1=colorbar; box on;
    caxis([minVal maxVal]);
    xlabel('x [cm]'); ylabel('y [cm]'); 
    title(['B_r at t = ',num2str(thist,3),' ns']);
    axis('equal'); 
    %axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);
    %
    hold on; 
    pf1=fill(xbox_ins, ybox_ins,[0.93 0.69 0.13]);
    hold off;
    axis([0 1.5 -0.1 1.2]);
    
    
    subplot(1,2,2);
    p1=pcolor(x,y,By_cc(:,:,it)); shading flat; colorbar; box on;
    xlabel('x [cm]'); ylabel('y [cm]'); 
    title(['B_\theta at t = ',num2str(thist,3),' ns']);
    axis('equal'); 
    %axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);
    %
    hold on; 
    pf1=fill(xbox_ins, ybox_ins,[0.93 0.69 0.13]);
    hold off;
    axis([0 1.5 -0.1 1.2]);
    %
% subplot(1,3,3);
% h23=pcolor(Yce,Xce,Ez(:,:,thisit)'); shading flat; colorbar; box on;
% xlabel('y [cm]');
% ylabel('x [cm]'); 
% title(['E_z at t = ',num2str(tns(thisit),3),' ns']);
% axis('equal'); 
% axis([Yce(2) Yce(end-1) Xce(2) Xce(end-1)]);


    map = colormap('jet');
    map(1,:) = [1 1 1]; % set lower bound to white
    colormap(map);


%     %%%   put time stamp on figure
%     %
%     a1=annotation(f1,'textbox',...
%     [0.0187 0.8775 0.085 0.08],...
%     'String',['t=',num2str(thist),'ns'],...
%     'FitBoxToText','off', ...
%     'backgroundcolor','y');


    % drawnow;
    %
    %
    %
    M(thisFigNum+1) = getframe(f1); 
    thisFigNum = thisFigNum+1;
    if(tns(it)~=tns(end))
       % delete(a1);
        delete(p1);
        delete(cb1);
        delete(pf1);
       % delete(p2);
       % delete(pf2);
    end


end


f1pos = get(f1,'position');
f11=figure(11);
set(f11,'Position',f1pos);
movie(f11,M,1,4);

v=VideoWriter('./liftOff2Drth.mp4', 'MPEG-4');
%v=VideoWriter('./dpf1Ddensity.mj2', 'Motion JPEG 2000');
%v=VideoWriter('./dpf1Ddensity.mj2', 'Archival');
%v=VideoWriter('./dpf1Ddensity.avi', 'Motion JPEG AVI');
%v=VideoWriter('./dpf1Ddensity.avi', 'UnCompressed AVI');
v.FrameRate = 4; v.Quality=100; %v.CompressionRatio = 2;
%v.LosslessCompression = true;
open(v); writeVideo(v,M);
close(v);
