clc
clear all

%declaration
num_points = 50;     %number of points per Rj (for disctasation of the integral and calc)

L0 = 12;
T0 = 6E6; %K;
beta = 2.6;
gamma = 2;
Fn = 330; %kg/s;
Fe = -2e6; %W; 
kb = 1.3806488E-23; %m^2*kg*s^-2*K^-1
Rj = 71492E3; %m
Z0 = 2;

%initalization of the Q map
[Q, L] = Qturbulence(num_points);
f = Q.*(L.^(1+gamma));
Integral_f_with_L = trapz(L,f);

L_length = length(L)

integral = zeros(1,L_length);
inside_integral=zeros(1, L_length);

%integration
for i=2:L_length  
    f_double_prime = f(1:i);
    L_double_prime = L(1:i);
    inside_integral(i) = trapz(L_double_prime,f_double_prime);
end

for i = 2:L_length
    L_prime_limit = i;
    
    L_prime = L(1:L_prime_limit);
    prime_integral=inside_integral(1:L_prime_limit);

    f_prime = L_prime.^(2-beta).*prime_integral;
    integral(i) = trapz(L_prime, f_prime);
end

%computation
test = (T0-Fe/(Fn))*(L/L0).^(beta-3-gamma);

% TA = (T0-Fe/(Fn))*(L/L0).^(beta-3-gamma)+Fe/(Fn)*(L/L0).^(-gamma);
TA = T0*(L/L0).^(-gamma);
TQ = 4*pi()*(3-beta)*Z0/(3*kb*Fn)*((L.^(beta-3-gamma)).*integral);

T = TA + TQ;

%ploting
figure(1)
semilogy(L,test)
title('test')
xlabel('Radial distance in Rj')
ylabel('T [K]')

figure(2)
semilogy(L,TQ)
title('turbulent heating')
xlabel('Radial distance in Rj')
ylabel('T [K]')

figure(3)
semilogy(L,TA)
title('cooling')
xlabel('Radial distance in Rj')
ylabel('T [K]')

figure(4)
semilogy(L,T)
title('Plasma temperature as a function of radial distance')
xlabel('Radial distance in Rj')
ylabel('T [K]')







