global Me gam Arat

Me = 0.9574;
N_Length = 10/39.37; %10 inches to meters
P300_Radius = .091/2;
x = linspace(0,N_Length);
y = P300_Radius-(x.^2/2);
A = pi.*y.^2;
P(1) = 134705.6523;
P0 = 179161.4235;
T = 1085.647611;
T0 = 1181.44485;
Arat=A(1)/A(end);
M(1) = LMFnlsq(@areaMach,0.3);
gamma(1)= 1.4 -8.32*10^-5 * (T(1)-273.15);

for g = [2:length(x)]
    Area_ratio = A(g)/A(end);
    gamma(g) = 1.4 -8.32*10^-5 * (T(g-1)-273.15);
    gam = gamma(g);
    Arat = Area_ratio;
    M(g) = LMFnlsq(@areaMach, M(1));
    P(g) = P0*(1+(gamma(g)-1)/2*M(g)^2)^(-gamma(g)/(gamma(g)-1));
    T(g) = T0 * (1+(gamma(g)-1)/2*M(g)^2)^-1;  
end

figure(1)
yyaxis left
plot (x*1000,P)
xlabel('position in nozzle (mm)') 
ylabel('Pressure (Pa)')

yyaxis right
plot (x*1000,M)
ylabel('Mach number') 

function [f] = areaMach(x)
    %Isentropic Flow Equations to solve for initial Mach given area ratio
    global Me gam Arat
    yf=gam;
    Mi=x(1);
    f=Me/Mi*sqrt(((1+(yf-1)/2*Me^2)/(1+(yf-1)/2*Mi^2))^((yf+1)/(yf-1)))-Arat;
end