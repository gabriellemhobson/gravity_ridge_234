% HW 7 Coding

% Gravity anomaly of a sediment covered ridge

clear all; close all;

n = 512; % grid refinement
Lx = 1e3; % length of x domain in m
Ly = 1e5; % length of y domain in m
Lz = 1e1; % length of z domain in m
x = linspace(-Lx/2,Lx/2,n);
y = linspace(-Ly/2,Ly/2,n);
z = linspace(0,Lz,n);

kx = (-n/2:n/2-1)/Lx;
ky = (-n/2:n/2-1)/Ly;
kz = (-n/2:n/2-1)/Lz;

sigma = 10*1000; % m
lambda_f1 = 5*1000; % flexural wavelength ?? (m)
lambda_f2 = 30*1000; % flexural wavelength ?? (m)
D1 = 1.8e21; % N m
D2 = 3.2e23; % N m
g = 9.82; 
G = 6.5e11; % gravitational constant, (m^3 kg^-1 s^-2)
d = 6*1000;
rho_w = 1025;
rho_s = 2300;
rho_c = 2800;
rho_m = 3330;
s = 6.2*1000;


% make the domain and wavenumber domain mesh
% [X,Y,Z] = meshgrid(x,y,z);
% [KX,KY,KZ] = meshgrid(kx,ky,kz);

% [X,Y,Z] = meshgrid(x,y,z);
[X,Y] = meshgrid(x,y);

% b = (2*1000)*exp(-(X.^2 + Y.^2 + Z.^2)/(2*sigma^2));  % where 2 km is the height of the ridge
b = (2*1000)*exp(-(X.^2)/(2*sigma^2));

figure
hold on
surf(X,Y,b)
shading flat

B = fftshift(fft2(fftshift(b)));

[KX,KY] = meshgrid(kx,ky);

k = sqrt(KX^2 + KY^2); % ??

R1 = 1./(1 + D1*abs(k).^4./(g*(rho_m-rho_c)));
R2 = 1./(1 + D2*abs(k).^4./(g*(rho_m-rho_s)));

Q2 = 2*pi*G*exp(-abs(k)*s) .* ( (rho_c - rho_s) + exp(-abs(k)*d).*(rho_m - rho_c).*(1+R2*(rho_s-rho_w)./(rho_m-rho_s)).^(-1).*(R2.*(rho_s-rho_w)/(rho_m-rho_s) - R1.*(rho_c-rho_w)./(rho_m-rho_c)));

GR = Q2.*B;

gr = ifftshift(ifft2(fftshift(GR)));

surf(X,Y,gr)
shading flat







