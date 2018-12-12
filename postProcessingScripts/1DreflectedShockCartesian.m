%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D dpf slug solution in cartesian
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

g = 5/3;   % adiabatic coefficient

c0 = g/2*(g-1);
c1 = (g^2-1)/4;
c3 = g-1;


a3 = 16*(c0-c1*(c3+1));
a2 = 16*c1*c3 - 8*(c0-c1*(c3+1));
a1 = c0 - c1*(c3+1)-8*c1*c3-1;
a0 = c1*c3;

x = roots([a3 a2 a1 a0]);

% %%%   plot piston and shock radius
% %
% close(figure(1));
% f1=figure(1); set(f1,'position',[850 11 560 780]);
% 




    