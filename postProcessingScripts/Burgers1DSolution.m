%%%%%%%%%%%%%
%%%
%%%   solution to Burger's equation from Hopf-Cole transformation
%%%
%%%   dF/dt + d(F^2/2)/dx = K*d2F/dx2
%%%
%%%%%%%%%%%%%
clear all;


%%%   set grid and initial profile
%
K = 0.002;
nx = 1000;
Xmin = -2.0;
Xmax = 2.0;
dx = (Xmax-Xmin)/nx;
x = Xmin-dx:dx:Xmax+dx;
a = 1;
b = 0.0;
c = 0.1;
xshift = (x-b)/c;
F0 = a*exp(-xshift.^2/2);


%%%   set time domain
%
time = 0:0.02:2.0;


%%%   initialize solution matrix and BC's
%
F = zeros(length(x),length(time));
F(:,1) = F0;
F(1,:) = F0(1);
F(end,:) = F0(end);


%%%   create integral vector of initial profile
%
arg2 = zeros(size(x));
[~,i0] = min(abs(x-0));
arg2(i0:end) = cumtrapz(x(i0:end),F0(i0:end));
for i=1:i0-1
    arg2(i) = -arg2(end+1-i);
end


%%%   loop over time
%
exparg = zeros(size(x));
lnFunc = zeros(size(F));
for n=2:length(time)
    for i=1:length(x)
        exparg = -(x(i)-x).^2/(4*K*time(n))-1/(2*K)*arg2;
        arg3=trapz(x,exp(exparg));
        lnFunc(i,n) = log(arg3/sqrt(4*pi*K*time(n)));
    end  
end


for n=2:length(time)
    for i=2:length(x)-1
        F(i,n) = -2*K*(lnFunc(i+1,n)-lnFunc(i-1,n))/2/dx;
    end  
end


%%%   plot solution
%
f1=figure(1);
plot(x,F(:,1));
hold on; plot(x,F(:,round(end/2)));
hold on; plot(x,F(:,end));
xlim([Xmin Xmax]);



