%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D dpf slug solution in cartesian
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


g=1.0:0.01:2.1;
xp0 = (g-1)./(g+1); % piston radius when shock reaches axis;
xs0 = (1-g).^2./g./(g+1); % piston/shock radius at stagnation
%
rhoshock = 1./xp0;     % ratio of shocked density to ambient
rhostag = 1./xs0;      % ratio of stag density to ambient
%
Ps0 = (3*g-1)./(g-1);     % ratio of stag pressure to mag pressure
%Ts0 = Ps0./(rhostag.^g);  % stag
%
usp = (g+1)/2.0;          % ratio of incident shock speed to piston speed
urs = 2*xp0;              % ratio of refl shock speed to inc shock speed


%%%   plot piston and shock radius
%
close(figure(11));
f11=figure(11);
set(f11,'position',[1000 845 560 500]);
%
subplot(2,2,1);
plot(g,xp0,'displayName','t_a');
hold on; plot(g,xs0,'displayName','t_s');
xlabel('\gamma'); ylabel('x_p/x_0');
grid on; axis([1 2 0 0.4]);
lg1=legend('show','location','northwest');
hold on; line([5/3 5/3],[0.0 0.25],'color','black','linestyle','--');
hold on; line([1 5/3],[0.25 0.25],'color','b','linestyle','--');
hold on; line([1 5/3],[0.1 0.1],'color','r','linestyle','--');
title('piston radii'); axis('square');
set(gca,'xtick',0:0.25:2);


%%%   plot ratio of refleted shock speed to incident shock speed
%%%   and incident shock speed wrt piston speed
%
subplot(2,2,3);
plot(g,usp,'displayName','U_i_s/U_p');
hold on; plot(g,urs,'displayName','U_r_s/U_i_s');
xlabel('\gamma'); %ylabel('U/U');
grid on; axis([1 2 0 1.6]);
lg2=legend('show','location','best');
hold on; line([5/3 5/3],[0.0 4/3],'color','black','linestyle','--');
hold on; line([1 5/3],[4/3 4/3],'color','b','linestyle','--');
hold on; line([1 5/3],[0.5 0.5],'color','r','linestyle','--');
title('speeds'); axis('square');
set(gca,'xtick',0:0.25:2);

    
%%%   plot ratio of stagnation pressure to piston pressure
%
subplot(2,2,4);
plot(g,Ps0);
xlabel('\gamma'); ylabel('P_s/P_B');
grid on; axis([1 2 0 20]);
hold on; line([5/3 5/3],[0.0 6],'color','black','linestyle','--');
hold on; line([1 5/3],[6 6],'color','b','linestyle','--');
title('stagnation pressure'); axis('square');
set(gca,'xtick',0:0.25:2);


%%%   plot densities
%
subplot(2,2,2);
plot(g,1./xp0,'displayName','post-shock');
hold on; plot(g,1./xs0,'displayName','stagnation');
xlabel('\gamma'); ylabel('\rho/\rho_0');
grid on; axis([1 2 0 20]);
lg4=legend('show','location','best');
hold on; line([5/3 5/3],[0.0 10],'color','black','linestyle','--');
hold on; line([1 5/3],[4 4],'color','b','linestyle','--');
hold on; line([1 5/3],[10 10],'color','r','linestyle','--');
title('densities'); axis('square');
set(gca,'xtick',0:0.25:2);



%%%   data from cylindrical simulations
%
g_cyl  = [1.2   1.4   5/3  2.0];
ra_cyl = [      0.005 0.28 0.40];
rs_cyl = [0.005 0.035 0.12 0.25];
%Ps_cyl = [          3.00];
figure(11); subplot(2,2,1);
hold on; plot(g_cyl,rs_cyl,'r*');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%   initial equation obtained by round about method. 
%%%   Two roots of cubic equation are always one.
%%%   general solutions are above

% g = 5/3;   % adiabatic coefficient
% 
% xp0 = (g-1)/(g+1); % piston radius when shock reaches axis;
% 
% c0 = g/2*(g-1);
% c1 = (g^2-1)/4;
% c3 = g-1;
% 
% a3 = 1/xp0^2*(c0-c1*(c3+1));
% a2 = 1/xp0^2*c1*c3 - 2/xp0*(c0-c1*(c3+1));
% a1 = c0 - c1*(c3+1)-2/xp0*c1*c3-1;
% a0 = c1*c3;
% 
% xs = roots([a3 a2 a1 a0]);  % rs/r0 at stagnation
% Ps = (g-1)*(1./xs-1);       % Ps/P0 at stagnation
% 
% display(xs');
