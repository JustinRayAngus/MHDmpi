function f0 = Poisson_Solver_2D(xplot0,zplot0,g0,plot_results,cyl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   Solve the 2D Poison equation with assumed period in z-direction
%%%   and zero value boundary conditions in the x-direction
%%%
%%%   (d2dx2 + d2dz2)f0(x,z) = g(x,z)
%%%   or
%%%   (d2dx2 + ddx/x + d2dz2)f0(x,z) = g(x,z)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5
    cyl = 0; % use cartesian by default
end

%display(cyl);

%%%   get rid of ghost cells
%
nx0 = length(xplot0(:,1));        % includes ghost cells
nz0 = length(zplot0(1,:));        % includes repeat cell
xplot = xplot0(3:end-2,3:end-2);  % cell center
zplot = zplot0(3:end-2,3:end-2);  % cell center
g     = g0(3:end-2,3:end-2);      % does not include ghost cells 


%%%   define grid stuff
%
nx = length(xplot(:,1));
nz = length(zplot(1,:));
dx = xplot(2,1)-xplot(1,1);
dz = zplot(1,2)-zplot(1,1);
Lz = 2*(max(max(zplot))+dz/2); % length of z-domain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   obtain fourier coefficients of g = sum gm(z)*cos(kmz)
%%%

if(rem(nz/2,1)==0)
    mmax = nz/2; % nz is even
else
    mmax = (nz-1)/2;
end
m  = 0:1:mmax-1;
km = 2*pi*m/Lz;
am = zeros(nx,mmax); % cosine components
bm = zeros(nx,mmax); % sine components

for i = 1:nx
    am(i,1) = sum(g(i,:))*dz/Lz;
    for n = 2:mmax
        am(i,n)  = 2*sum(g(i,:).*cos(km(n)*zplot(i,:)))*dz/Lz;
        bm(i,n) = 2*sum(g(i,:).*sin(km(n)*zplot(i,:)))*dz/Lz;
    end
end


gtest = zeros(size(g));
for j = 1:length(zplot(1,:))
    gtest(:,j) = am*cos(km'*zplot(1,j));
    gtest(:,j) = gtest(:,j) + bm*sin(km'*zplot(1,j));
end
%figure(11);
%error = abs(g-gtest);
%[~,h] = contourf(zplot,xplot,error,30); axis('equal'); colorbar;
%[~,h] = contourf(zplot,xplot,gtest,30); axis('equal'); colorbar;
%set(h,'linestyle','none');
%plot(zplot(1,:),g(130,:));
%hold on; plot(zplot(1,:),gtest(130,:),'r--');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%       Thomas algorithm  for (d2dx2 - km^2)fm(x) = gm(x)
%%%

lbc = -1;   % lower bc, -1 = odd, 0 = 0 at ghost cell, 1 = even
ubc = 0;    % upper bc, -1 = odd, 0 = 0 at ghost cell, 1 = even
fa = zeros(size(am)); % cosine coefficients of solution
fb = zeros(size(bm)); % sine coefficients of solution

for n = 1:mmax
    
    a = ones(nx,1);
    c = ones(nx,1);
    if(cyl)
        a = a - dx/2.0./xplot(:,1);
        c = c + dx/2.0./xplot(:,1);
    end
    b = -(2+km(n)^2*dx^2)*ones(nx,1);  % diagonal compoments
    b(1)  = b(1) + lbc; % upper bc
    b(nx) = b(nx)+ ubc; % lower bc
    ha = am(:,n)*dx^2;  % A*f = h
    hb = bm(:,n)*dx^2;
   
    %%%    reduce to upper tridiagonal matrix
    %
    for i=2:nx
        b(i)  = b(i)  - a(i)/b(i-1)*c(i-1);
        ha(i) = ha(i) - a(i)/b(i-1)*ha(i-1);
        hb(i) = hb(i) - a(i)/b(i-1)*hb(i-1);
    end
    
    %%%    back solve
    %
    fa(nx,n) = ha(nx)/b(nx);
    fb(nx,n) = hb(nx)/b(nx);
    for i=nx-1:-1:1
        fa(i,n) = (ha(i)-c(i)*fa(i+1,n))/b(i);
        fb(i,n) = (hb(i)-c(i)*fb(i+1,n))/b(i);
    end
    
end

%%%     create solution by summing over fourier compomnents
%
f=zeros(size(g));
for j = 1:length(zplot(1,:))
    f(:,j) = fa*cos(km*zplot(1,j))' + fb*sin(km*zplot(1,j))';
end
% [~,h] = contourf(zplot,xplot,f,30); axis('equal'); colorbar;
% set(h,'linestyle','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%      compute Laplacian to see if solution is correct
%%%

%%%    extend solution to domain with x and z ghost cells
%
f0 = zeros(size(xplot0));
f0(3:nx0-2,3:nz0-2) = f;
f0(:,2) = f0(:,end-2);
f0(:,1) = f0(:,end-3);
f0(:,end) = f0(:,4);
f0(:,end-1) = f0(:,3);
%
f0(2,:) = lbc*f0(3,:);
f0(1,:) = lbc*f0(4,:);
f0(nx0-1,:) = ubc*f0(nx0-2,:);
f0(nx0,:) = ubc*f0(nx0-3,:);


%%%    compute d2dz2f (periodic)
%
d2fdz2 = zeros(size(f));
d2fdz2(:,1) = (f(:,nz)-2*f(:,1)+f(:,2))/dz^2;
for j = 2:nz-1
    d2fdz2(:,j) = (f(:,j-1)-2*f(:,j)+f(:,j+1))/dz^2;
end
d2fdz2(:,nz) = (f(:,nz-1)-2*f(:,nz)+f(:,1))/dz^2;


%%%    compute d2fdx2 and dfdx
%
d2fdx2 = zeros(size(f));
dfdx = zeros(size(f));
for i = 1:length(xplot)-1
    d2fdx2(i,:) = (f0(i-1+2,1:nz)-2*f0(i+2,1:nz)+f0(i+1+2,1:nz))/dx^2;
    dfdx(i,:)   = (f0(i+1+2,1:nz)-f0(i-1+2,1:nz))/2.0/dx;
end

%%%  now add together to create laplacian
%
g2 = d2fdx2 + d2fdz2;
if(cyl)
    g2 = g2 + dfdx./xplot;
end
%errorg = abs(g-g2);

%figure(7); plot(zplot(130,:),g2(20,:)); grid on;

if(plot_results==1)
    %%% plot both g, g2, and error
    %
    close(figure(3));
    f3=figure(3);
    scrsz = get(0,'ScreenSize');
    set(f3,'Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2]);
    %
    subplot(1,3,1); 
    [~,ha] = contourf(zplot,xplot,g,30); axis('equal'); colorbar;
    set(ha,'linestyle','none');
    xlabel('z/\delta'); ylabel('x/\delta'); title('Original g');
    %
    subplot(1,3,3); 
    [~,ha] = contourf(zplot,xplot,g2,30); axis('equal'); colorbar;
    set(ha,'linestyle','none');
    xlabel('z/\delta'); ylabel('x/\delta'); title('g from laplacian');
    %
    subplot(1,3,2); 
    [~,ha] = contourf(zplot,xplot,f,30); axis('equal'); colorbar;
    set(ha,'linestyle','none');
    xlabel('z/\delta'); ylabel('x/\delta'); title('f from \nabla^2f=g');
    %
end


end