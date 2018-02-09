%%%%%%%%%%%%%%%%%%%%
%%%
%%%    d(P+B^2/2)/dr + B^2/r  = 0
%%%    or 
%%%    d(rP)/r - P/r + d(rB)/dr/r*B = 0;
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
%P(end) = 0.010185410965725;
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
g = rcc.*P;
h = rcc.*B;
dPdr = zeros(size(rcc));
dhdr = zeros(size(rcc));
Jcc =  zeros(size(rcc));

for i=2:length(rcc)-1
    dPdr(i) = (P(i+1)-P(i-1))/2/dr;
    Jcc(i)  = (h(i+1)-h(i-1))/2/dr/rcc(i);  
end
%P(end) = P(end-2)-Jcc(end-1)*B(end-1)*2*dr;

Fr0 = -(dPdr + Jcc.*B);
Fr0(1) = 0;
Fr0(end) = 0;

figure(3); hold on; plot(rcc,Fr0); box on; grid on;
title('radial force');

%hold on; plot(rcc,Fr0./P);
%figure(3); hold on; plot(rcc,Fr1,'rx');
%figure(2); hold on; plot(rcc,dfdr);
%hold on; plot(rcc,dB2over2dr,'r--');
%hold on; plot(rcc,dPdr,'r--');
%hold on; plot(rcc,dPdr+dB2over2dr,'r--');



%%%   calculate Pressure s.t force balance exactly at cell-edges
%
Jce = zeros(size(rce));
Bce = zeros(size(rce));
Ptest = zeros(size(rcc));
Ptest(1:2) = P(1:2);

for i=2:length(rcc)
   Jce(i-1) = (rcc(i)*B(i) - rcc(i-1)*B(i-1))/dr/rce(i-1);
   Bce(i-1) = (B(i)+B(i-1))/2.0;
end
Jce(1) = sqrt(2)*2/a;

for i=3:length(rcc)-1
    Ptest(i) = Ptest(i-1) - Jce(i-1)*Bce(i-1)*dr;
end
Ptest(end)=Ptest(end-1);


figure(1); hold on; plot(rcc,Ptest,'g--');


%%%   recalculate force balance using Ptest
%
dPtestdr = zeros(size(rcc));
dPtestdrce = zeros(size(rce));
Jcc = zeros(size(rcc));
for i=2:length(rcc)-1
    dPtestdr(i) = (Ptest(i+1)-Ptest(i-1))/2/dr;  
    dPtestdrce(i) = (Ptest(i+1)-Ptest(i))/dr;  
    Jcc(i) = (Jce(i)+Jce(i-1))/2.0;
end
Jcc(end) = Jce(end);
%dPtestdrce(1) = 0;

FrTestce = dPtestdrce + Jce.*Bce;
%FrTest = dPtestdr + Jcc.*B;
FrTest = dPdr + Jcc.*B;

%figure(3); hold on; plot(rce,FrTestce,'g');
%figure(3); hold on; plot(rcc,FrTest,'r--');


