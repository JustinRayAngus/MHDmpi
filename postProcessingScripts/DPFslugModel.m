%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D dpf slug solution from potter 1978
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

gamma = 5/3;

dt = 0.0001;
tmax = 0.4;
t = 0:dt:tmax;

Vs = zeros(size(t));   % shock speed
Vp = zeros(size(t));   % piston speed
rs = zeros(size(t));   % shock radius
rp = zeros(size(t));   % piston radius


%%%   set initial condition
%
rs(1) = 1;
rp(1) = 1;
Vs(1) = -sqrt(1+gamma)/rp(1);


%%%   march in time until rs=0
%
for i=2:length(t)
   
    %%%   predictor stage
    %
    thisdt = dt/2;
    
    drs = Vs(i-1)*thisdt;
    factp = (gamma-1)/gamma + rs(i-1)^2/rp(i-1)^2/gamma;
    facts = 2/(gamma+1)*rs(i-1)/rp(i-1);
    drp = facts/factp*drs;
    
    rs_half = rs(i-1) + drs;
    rp_half = rp(i-1) + drp;
    Vs_half = -sqrt(1+gamma)/rp_half;
    
    %%%   corrector stage
    %
    thisdt = dt;
    drs = Vs_half*dt;
    factp = (gamma-1)/gamma + rs_half^2/rp_half^2/gamma;
    facts = 2/(gamma+1)*rs_half/rp_half;
    drp = facts/factp*drs;
    
    rs(i) = rs(i-1) + drs;
    rp(i) = rp(i-1) + drp;
    Vs(i) = -sqrt(1+gamma)/rp(i);
    
    if(rs(i)<=0)
        rs(i+1:end) = zeros;
        rp(i+1:end) = rp(i);
        Vs(i+1:end) = Vs(i);
        itmax = i;
        break;
    end

end
x = rs./rp;
factp = (gamma-1)/gamma + x.^2/gamma;
facts = 2/(gamma+1)*rs./rp;
Vp = facts./factp.*Vs; % piston velocity
Vfs = 2/(gamma+1)*Vs; % post-shock fluid velocity


%%%   set analytical solution for rp from Eq. 12
%
rp_soln = (gamma./(gamma+1-x.^2)).^(gamma/(gamma-1));
ap = (gamma/(gamma+1))^(gamma/(gamma-1)); % final piston radius


%%%   plot piston and shock radius
%
close(figure(1));
f1=figure(1); set(f1,'position',[850 11 560 780]);

subplot(2,1,1);
plot(t,rp,t(1:itmax),rs(1:itmax)); title('positions');
hold on; line([t(1) t(itmax)],[ap ap],'linestyle','--','color','black');
xlabel('t/t_0'); ylabel('r/r_0');
lg1=legend('piston','shock'); set(lg1,'location','best');
grid on; axis([t(1) t(end) 0 1]);


subplot(2,1,2);
plot(t,Vp,t,Vs); title('velocities');
hold on; plot(t,Vfs,'linestyle','--');
xlabel('t/t_0'); ylabel('V/V_0');
lg1=legend('piston','shock','post-shock fluid'); set(lg1,'location','best');
grid on; %axis([t(1) t(end) 0 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   calculate stuff from expressions used by Potter
%%% 


%%%   calculate volume and pressure
%
Volume = pi*(rp.^2-rs.^2);
Pressure = 2/(gamma+1)*Vs.^2;


%%%   calculate dVol/dt used by Potter where post shock fluid speed
%%%   is used for shock front location
%
dVol_dt_potter = 2*pi*(rp.*Vp - rs.*Vfs); 
dVol_dt_actual = 2*pi*(rp.*Vp - rs.*Vs); 


f3=figure(3); 
plot(t,dVol_dt_potter,'displayName','adiabatic (potter)');
hold on; plot(t,dVol_dt_actual,'displayName','total');
title('dVolume/dt'); xlabel('t/t_0'); ylabel('dVol/dt');
grid on;
lg3=legend('show','location','best');



    