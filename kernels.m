%%%% Compute the different bounds for acceptable kernel examples, in the paper "A Dual Certificates Analysis of Compressive Off-the-Grid Recovery"
%%%% contact author: nicolas.keriven@gmail.com

clear all

%% for Gaussian kernel, with separation sqrt(2*(a*log(k)+b*log(d)+c))*sig

a = 5;
b = 2;
c = 12;
e = exp(1);

deltaGauss = exp(-c)*(1+4*a/(e*(a-1)) + 4*b/(e*(b-1/2)) + 4*c)

delta2Gauss = exp(-c) + 2*exp(-2*c)/(1-deltaGauss)*(a/(2*e*(a-1)) + b/(e*(2*b-1)) + c)

epsGauss = 1-((0.7789+exp(-c/4))/(1-delta2Gauss) + sqrt(2)*exp(-c)*sqrt(a/(2*e*(a-1))+ b/(e*(2*b-1)) + c)/((1-deltaGauss)*(1-delta2Gauss))*(0.6066 +...
    exp(-c/4)/sqrt(2)*sqrt(a/(e*(a/2-2)) + 2/e + c)))

lambdaGauss = (1 - delta2Gauss/(1-delta2Gauss))*0.3893 - (1 + a/(e*(a-4)) + 1/e + c/4)*exp(-c/4)/(1-delta2Gauss) -...
    sqrt(2)*exp(-c)*sqrt(a/(2*e*(a-1))+ b/(e*(2*b-1)) + c)/((1-deltaGauss)*(1-delta2Gauss))*(2.4750+...
    exp(-c/4)/(2*sqrt(2))*(6*(a/(e*(a/2-2)) + 2/e + c)^(1/2) + (3*a/(e*(a/2-2)) + 6/e + c)^(3/2)))


%% multi-dimensional Fejer kernel

%%% bound for t < 0.7559 lambda_c
eps0 = @(c)(pi^2*c^2/6 - pi^4*(1+1/64)^4*c^4/72); % times 1/d, for t = c/(fc*sqrt(d)). K0near(t) = 1 - eps0
K1near = @(c)(pi^2*(1+1/32)*c/3); % times fc, t = c/fc
K2near = pi^3*(1+1/32)/3; %  times fc^2
K3near = @(c)(pi^4*(1+1/64)^4*c/3);
epslower = @(c)(pi^2*c^2/6);

%%% bound for far
%%% Uniform bounds for t = C*d^(1/2)*k^(1/4)/fc, for fc >
%%% 128*d^(1/2)*k^(1/4) and C > 0.7559
a = @(C)(2/(pi*(1-pi^2*C^2/(6*128^2))));
b = @(C)(a(C)/C);
K0far = @(C)(a(C)^4 / C^4); % times 1/(kd^2)
K1far = @(C)(pi*a(C)^4*(2 + 2*b(C))/C^4); % times fc/(kd^2)
K2far = @(C)(pi^2*a(C)^4*(4 + 7*b(C) + 6*b(C)^2)/C^4); % times fc^2/(kd^2)
K3far = @(C)(pi^3*a(C)^4*(8+24*b(C)+30*b(C)^2+15*b(C)^3)/C^4); % times fc^3/(kd^2)

K1max = max(K1near(0.7559),K1far(0.7559));
K3max = max(K3near(0.7559),K3far(0.7559));
K2max = K2near;

%%% parameters
Cnear = 0.1; % taking eps_near = Cnear / (sqrt(d)*fc)
Delt = 5; % taking Delta = Delt*d^(1/2)*k^(1/4) / fc
%%% intermediate bounds
v = pi^2/3;
lambda1 = (pi^2/3 - pi^4*(1+1/64)^4*Cnear^2)*exp(-(epslower(Cnear)+1)/(1-epslower(Cnear))) - K1near(Cnear)^2; % times fc^2. Because (1-eps/d)^d <= e^(-eps)
b3 = K3near(Cnear) + K1near(Cnear)*K2near + K1near(Cnear)^3; % times d*fc^3
e2 = K2far(Delt/2) + K2max*K0far(Delt/2) + K1max*K1far(Delt/2)+K1max^2*K0far(Delt/2); % times fc^2/(kd)
e3 = K3far(Delt/2) + K3max*K0far(Delt/2) + K1max*K2max*K0far(Delt/2) + K2max*K1far(Delt/2) + K1max*K2far(Delt/2); % times fc^3*sqrt(d)
%%% final bounds
u = K1far(Delt)/sqrt(v); % times 1/d^2
deltaFejer = K2far(Delt) % times 1/d
delta2Fejer = K0far(Delt) + u^2*(1-deltaFejer)^(-1) % times 1/d^2
epsFejer = (eps0(Cnear)-K0far(Delt/2))*(1-delta2Fejer)^(-1) - delta2Fejer/(1-delta2Fejer) - K1far(Delt)*(K1max + K1far(Delt/2))/(sqrt(v)*(1-deltaFejer)*(1-delta2Fejer)) % times 1/d
lambdaFejer = (1-delta2Fejer/(1-delta2Fejer))*lambda1 - e2/(1-delta2Fejer) - K1far(Delt)*(b3 + e3)/(sqrt(v)*(1-deltaFejer)*(1-delta2Fejer)) % times fc^2
