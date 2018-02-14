%%%%%%%%%%%%%%%%%%%%
%%%
%%%    d(P+B^2/2)/dr + B^2/r  = 0
%%%    or 
%%%    d(rP)/r - P/r + d(rB)/dr/r*B = 0;
%%%
%%%    d(r*(P+B^2/2))/dr - (P-B^2/2) = 0
%%%
%%%%%%%%%%%%%%%%%%%%
clear all;

nr = 200;
R  = 1;
a  = R/3;
dr = R/nr;

rcc = -dr/2:dr:R+dr/2;
rce = rcc(1:end-1)+dr/2;
%rcc = abs(rcc);

x = rcc/a;
P = 1./(1+x.^2).^2;
B = sqrt(2)*x./(1+x.^2);
P(end) = P(end-1);
B(end) = B(end-1)*rcc(end-1)/rcc(end);
%P(end) = 2*P(end-1) - P(end-2);
%B(end) = 2*B(end-1) - B(end-2);
%P(end) = 3*(P(end-1) - P(end-2)) + P(end-3);
%B(end) = 3*(B(end-1) - B(end-2)) + B(end-3);

dPdr = -4*x./(1+x.^2).^3/a;
dB2over2dr = (-4*x.^3./(1+x.^2).^3+2*x./(1+x.^2).^2)/a;

figure(1); plot(rcc,P,'black');
hold on; plot(rcc,B,'blue');


%%%   calculate numerical force
%
f = P+B.^2/2;
Fluxcc = rcc.*f;
divFlux = zeros(size(rcc));

for i=2:length(rcc)-1
    divFlux(i) = (Fluxcc(i+1)-Fluxcc(i-1))/2.0/dr;
end

Force = -(divFlux - f + B.^2)./rcc;
Force(1) = 0;
Force(end) = 0;
Force(end-1) = Force(end-2)/3.0;
figure(3); plot(rcc,Force);


%%%   redefine magnetic field from force balance
%
B = sqrt(f-divFlux);
B(1) = -B(2);
B(end) = B(end-1)*rcc(end-1)/rcc(end);

figure(1); hold on; plot(rcc,B,'r--');


%%%   redefine pressure 
%
P = f-B.^2/2;
P(1) = P(2);
P(end) = P(end-1);

figure(1); hold on; plot(rcc,P,'r--');


%%%   recalculate numerical force
%
f = P+B.^2/2;
Fluxcc = rcc.*f;

for i=2:length(rcc)-1
    divFlux(i) = (Fluxcc(i+1)-Fluxcc(i-1))/2.0/dr;
end

Force = -(divFlux - f + B.^2)./rcc;
Force(1) = 0;
Force(end) = 0;
Force(end-1) = Force(end-2)/3.0;
figure(3); hold on; plot(rcc,Force);
legend('numerical force using C2','numerical force with redefined B and P');


