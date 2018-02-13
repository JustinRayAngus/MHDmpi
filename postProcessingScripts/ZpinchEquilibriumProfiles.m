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

nr = 100;
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

figure(1); plot(rcc,P,'b');
hold on; plot(rcc,B,'r');


%%%   calculate numerical force balance using C2
%
f = P;
g = rcc.*(P+B.^2/2);
divFlux = zeros(size(rcc));
h = rcc.*B;
Jcc =  zeros(size(rcc));
Pcc =  zeros(size(rcc));
Bcc =  zeros(size(rcc));

for i=2:length(rcc)-1
    divFlux(i) = (g(i+1)-g(i-1))/2/rcc(i)/dr;
    Jcc(i)  = (h(i+1)-h(i-1))/2/dr/rcc(i);  
    Pcc(i)  = (P(i+1)+P(i-1))/2.0;
    Bcc(i)  = (B(i+1)+B(i-1))/2.0;
end
Pcc(1) = Pcc(2);
Pcc(end) = Pcc(end-1);
Pcc(end) = 2*Pcc(end-1)-Pcc(end-2);

Fr0 = -divFlux + (P-B.^2/2)./rcc;
Fr0(end-1) = Fr0(end-2)/3;
Fr0(1) = 0;
Fr0(end) = 0;

figure(3); hold on; plot(rcc,Fr0); box on; grid on;
title('radial force');


%%%   calculate pressure such that force balance is exact assuming C2
%
Ptest = zeros(size(rcc));
Ptest(1:2) = P(1:2);
%Ptest(1:2) = Pcc(1:2); % works better at r=dr/2 for some resolutions?

for i=2:length(rcc)-1
    Ptest(i+1) = -rcc(i+1)*B(i+1)^2/2 + rcc(i-1)*(Ptest(i-1) + B(i-1)^2/2) + 2*dr*(Ptest(i)-B(i)^2/2);
    Ptest(i+1) = Ptest(i+1)/rcc(i+1);
end
Ptest(end) = Ptest(end-1);
    
figure(1); hold on; plot(rcc,Ptest,'r--'); grid on;


%%%   recalculate force
%
g1 = rcc.*(Ptest+B.^2/2);
divFlux1 = zeros(size(rcc));


for i=2:length(rcc)-1
    divFlux1(i) = (g1(i+1)-g1(i-1))/2/rcc(i)/dr;
end

Fr1 = -divFlux1 + (Ptest-B.^2/2)./rcc;
Fr1(end-1) = Fr1(end-2)/3;
Fr1(1) = 0;
Fr1(end) = 0;

figure(3); hold on; plot(rcc,Fr1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   calculate pressure from force balance using Newton itertions
%%%   (this will be useful when using TVD for flux rather than C2)
%%%
%%%

f = P+B.^2/2;
Fluxcc = rcc.*f;

for i=2:length(rcc)-1
    divFlux(i) = (Fluxcc(i+1)-Fluxcc(i-1))/2.0/dr;
end
divFlux(1) = divFlux(2);

Force = (divFlux - f + B.^2)./rcc;
Force(1) = 0;
Force(end) = 0;
Force(end-1) = Force(end-2)/3.0;
%divFlux = rcc.*Force+f-B.^2;
figure(3); hold on; plot(rcc,Force,'g');


%%%   redefine magnetic field and recalculate force
%
B = sqrt(f-divFlux);
B(1) = -B(2);
B(end) = B(end-1)*rcc(end-1)/rcc(end);

% f = P+B.^2/2;
% Fluxcc = rcc.*f;
% for i=2:length(rcc)-1
%     divFlux(i) = (Fluxcc(i+1)-Fluxcc(i-1))/2.0/dr;
% end
% 
% Force = (divFlux - f + B.^2)./rcc;
% Force(1) = 0;
% Force(end) = 0;
% Force(end-1) = Force(end-2)/3.0;
% 
% figure(3); hold on; plot(rcc,Force,'black');



