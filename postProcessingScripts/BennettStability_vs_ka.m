%%%%%%%%%%%%%%%%%%%%
%%%
%%%    d(P+B^2/2)/dr + B^2/r  = 0
%%%    or 
%%%    d(rP)/r - P/r + d(rB)/dr/r*B = 0;
%%%
%%%    d(r*(P+B^2/2))/dr - (P-B^2/2) = 0
%%%
%%%    calculate stability region from Kadomtsev 1960 Eqs 29 and 31
%%%    for Bennett equilibrium profiles of zpinch
%%%
%%%    P = a^4/(a^2+r^2)^2
%%%    B = sqrt(2)*a*r/(a^2+r^2)
%%%
%%%    gamma = 5/3, beta = 2*P/B^2
%%%
%%%    ideal sausage: w^2 ~ dlnP/dlnr + 4*gamma/(2+gamma*beta)
%%%
%%%    drift-ideal entropy (if sausage satisifed)
%%%    dlnT/dlnr + 1 + (7/10 + beta/4)*dlnP/dlnr > 0:
%%%    Note that 7/10 = (2*gamma-1)/(2*gamma), gamma = 5/3
%%%
%%%
%%%    r is normalized by some r0, 
%%%    U is normalized by sqrt(P0/(MiN0)) = sqrt(2*T0/Mi) = Cs0
%%%    omega is normalized by omega0=Cs0/r0
%%%
%%%%%%%%%%%%%%%%%%%%
clear all;


%%%   set radial grid
%
nr = 400;
R  = 3;
a  = 1;
dr = R/nr;


%%%   calculate Bennett profiles, qP, and plasma beta
%
rcc = dr/2:dr:R-dr/2;
x = rcc/a;
P = 1./(1+x.^2).^2;
B = sqrt(2)*x./(1+x.^2);
P(end) = P(end-1);
B(end) = B(end-1)*rcc(end-1)/rcc(end);
qP = -4*x.^2./(1+x.^2); % dlnP/dlnr
qT = 0*x;               % dlnT/dlnr
beta = 1./x.^2;         % plasma beta, 2*P/B^2;



%%%   specify input parameters
%
ka0 = 1:1:30;
growthRate_entropy0 = zeros(size(ka0));
growthRate_ideal0   = zeros(size(ka0));
growthRate_general0 = zeros(length(ka0),4);
gamma = 5/3; %2.1; %5/3; %2.1;  % adiabatic coefficient
%ka = 3.0;      % normalized wavenumber
%Li = 1.5e-2;  % omega0/Omegai0 = ion inertial length / r0
Li = 0.12;
nuT = 100; %100;    % normalized thermalization rate



%growthRate_ideal0 = zeros(size(ka0));
%growthRate_entropy0 = zeros(size(ka0));
%growthRate_ext0 = zeros(size(ka0));
%growthRate_total0 = zeros(size(ka0));
for k=1:length(ka0)
    ka = ka0(k);
    krhoi = ka*Li./B/sqrt(2);
    krhos = krhoi*sqrt(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%%
    %%%          compute modes in limit where k*rhos = kR*Li << 1
    %%%          Kadomtsev 1960 Eqs 29 and 30
    %%%
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%   compute growth rate from ideal MHD (neglect magnetosonic waves)
    %
    omega2 = 2.0./rcc.^2.*(qP + 4*gamma./(2 + gamma*beta));
    omega = sqrt(omega2);

    growthRate_ideal0(k)   = max(max(imag(omega)),0);


    %%%   compute growth rate from entropy mode (Kadomtsev)
    %
    nuTstar = nuT*rcc.*B/(ka*Li);
    a0 = (2*rcc.*B/(ka*Li)).^2;
    b0 = -2*rcc.*B/(ka*Li).*(2*qP - 1i*nuTstar);
    RHS = gamma^2*(4 + beta.*qP).^2./(4*gamma + qP.*(2 + gamma*beta)).* ...
          (1 + ((2*gamma-1)./(2*gamma)+beta/4).*qP - qT);
    c0 = qP.*(qP -4i*nuTstar) - RHS;

    omega_entropy = zeros(length(rcc),2);
    for i=1:length(rcc)
        thisp = [a0(i) b0(i) c0(i)];
        omega_entropy(i,:) = roots(thisp);
    end

    growthRate_entropy0(k) = max( max(imag(omega_entropy(:,1))), ...
                                  max(imag(omega_entropy(:,2))) );
    growthRate_entropy0(k) = max(growthRate_entropy0(k),0);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%%
    %%%          compute modes from MHD with full Ohms law, but
    %%%          only standard entropy density equation
    %%%
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    c3_ext = (gamma-1)/gamma - (1+beta/2);
    c2_ext = -(gamma-1)/gamma*qP - 4 + qT;
    c1_ext = 4./krhoi.^2.*(2 + (1+beta/2).*qP - (gamma-1)/gamma*qP);
    c0_ext = 4*qP./krhoi.^2.*((gamma-1)/gamma*qP - qT);

    omega_ext = zeros(length(rcc),3);
    for i=1:length(rcc)
        thisp = [c3_ext(i) c2_ext(i) c1_ext(i) c0_ext(i)];
        omega_ext(i,:) = roots(thisp);
        omega_ext(i,:) = omega_ext(i,:)*ka*Li/(2*rcc(i)*B(i));
    end
    
    
    growthRate_ext0(k) = max( max(imag(omega_ext(:,1))), ...
                              max(imag(omega_ext(:,2))) );
    growthRate_ext0(k) = max(growthRate_ext0(k),max(imag(omega_ext(:,3))));
    growthRate_ext0(k) = max(growthRate_ext0(k),0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%%
    %%%          compute modes from full dispersion (still no magnetosonic)
    %%%          Kadomtsev 1960 Eqs 26-28 (with corrections)
    %%%
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%   this is completely general relation with thermalization and qT
    %
    a = 2 + ((gamma-1)/gamma + beta/2).*qP - qT;
    b = 8./krhos.^2.*((gamma-1)/gamma.*qP - qT);
    c = gamma*(2+beta/2.*qP);
    d = gamma*beta/2.*qT;
    e = qP-qT;
    f = 8./krhos.^2.*(2+(1+beta/2).*qP-qT);
    g = (1+beta/2);
    h = qT-4;
    c4new = a.*(g - (gamma-1)/gamma);
    c3new = e.*qT-a.*h + (a.*g-e*(gamma-1)/gamma).*(4i.*nuTstar - 2*qP) ...
          - qP.*(e-a)*(gamma-1)/gamma - (e-a).*((gamma-1)/gamma*(qP-4i.*nuTstar)+qT);
    c2new = e.*b-a.*f + (4i.*nuTstar - 2*qP).*(e.*qT-a.*h) ...
          - ((4i.*nuTstar-qP).*qP + a.*c).*(a.*g - e*(gamma-1)/gamma) ...
          - (e-a).*(a.*d-qP.*qT+4i.*nuTstar.*qT+b) ...
          + (e-a).*qP.*((gamma-1)/gamma*(qP-4i*nuTstar)+qT);
    c1new = (4i.*nuTstar-2*qP).*(e.*b-a.*f) - ((4i.*nuTstar-qP).*qP+a.*c).*(e.*qT-a.*h) ...
          - (e-a).*b.*(4i.*nuTstar-qP) ...
          + qP.*(e-a).*(a.*d-qP.*qT+4i*nuTstar.*qT+b);
    c0new = -(qP.*(4i.*nuTstar-qP)+a.*c).*(e.*b-a.*f) + qP.*(e-a).*b.*(4i.*nuTstar - qP);
    
    
    %%%   expressions in limit of no thermalization and no qT
    %
    c4 = (gamma-1)/gamma-(1+beta/2);
    c3 = 2*qP.*(1+beta/2-(gamma-1)/gamma) - 4;
    c2 = (gamma-1)/gamma*qP.*(qP-gamma*(2+beta/2.*qP)-4./krhoi.^2) ... 
       + 4./krhoi.^2.*(2+(1+beta/2).*qP) + 8*qP ...
       - (1+beta/2).*(qP.^2-gamma*(2+beta/2.*qP).*(2+((gamma-1)/gamma+beta/2).*qP));
    c1 = (gamma-1)/gamma*8*qP.^2./(krhoi).^2 - 8*qP.*(2+(1+beta/2).*qP)./krhoi.^2 ...
       - 4*qP.^2 + 4*gamma*(2+beta/2.*qP).*(2+((gamma-1)/gamma+beta/2).*qP);
    c0 = -4*qP.^2./krhoi.^2.*(gamma-1)/gamma.*(qP-gamma*(2+beta/2.*qP)) ...
       + 4./krhoi.^2.*(2+(1+beta/2).*qP).*(qP.^2-gamma*(2+beta/2.*qP).*(2+((gamma-1)/gamma+beta/2).*qP));
    
    % c4 = c4.*(2*rcc.*B/(kR*Li)).^4;
    % c3 = c3.*(2*rcc.*B/(kR*Li)).^3;
    % c2 = c2.*(2*rcc.*B/(kR*Li)).^2;
    % c1 = c1.*(2*rcc.*B/(kR*Li)).^1;
    
    
    omega_total = zeros(length(rcc),4);
    for i=1:length(rcc)
     %   thisp = [c4(i) c3(i) c2(i) c1(i) c0(i)];
        thisp = [c4new(i) c3new(i) c2new(i) c1new(i) c0new(i)];
        omega_total(i,:) = roots(thisp);
        omega_total(i,:) = omega_total(i,:)*ka*Li/(2*rcc(i)*B(i));
    end


    for m=1:4
       % [~,ir0]=min(abs(rcc-1)); % where r=a
       % growthRate_general0(k,m) = max(imag(omega_total(ir0,m)),0);
        growthRate_general0(k,m) = max( max(imag(omega_total(:,m))),0);
    end
    growthRate_total0(k) = max( growthRate_general0(k,:) );

    
end

%close(figure(1));
f1=figure(1);
hold on; plot(ka0,growthRate_ideal0,'linestyle','--','displayName','ideal');
hold on; plot(ka0,growthRate_entropy0,'linestyle','--','displayName','entropy');
hold on; plot(ka0,growthRate_ext0);
hold on; plot(ka0,growthRate_total0,'displayName','general');
title('maximum growth rate'); box on;
axis('square'); grid on; xlabel('ka'); 
ylabel('growth rate ');
%axis([min(ka0) max(ka0) 0 1.2]);
