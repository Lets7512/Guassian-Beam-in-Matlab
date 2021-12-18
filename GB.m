%syms x y z;
%Inputs
lambda = 3.4e-6; % defined per group (Group 8)
wo = 30e-6; % defined per section (section 3)
k = 2*pi/lambda;
A = 1;
zo = pi*(wo^2)/lambda;
step = (2^0.5)*pi/k;
x=[-5*wo:step:-step 0 step:step:5*wo];
y=[-5*wo:step:-step 0 step:step:5*wo];
%z=linspace(-zo*1.5,zo*1.5,1000);
%x = -wo*5:step:wo*5;
%y = -wo*5:step:wo*5;
%Calculations
L= length(x);
%k = [kx;ky;kz];
%r = [x;y;z];
%E = @(var) real(A.*exp(-1i.*var));
p = @(x,y) (x.^2+y.^2);
U = @(x,y) A.*exp(-p(x,y)/(wo^2));
%z = 0*wo;
[X,Y] = meshgrid(x,y);
%------------------------------------------
%Part 1
%------------------------------------------
figure(1)
subplot(2,1,1);
GB=U(X,Y).^2;
imagesc(x,y,GB)
caxis([0 1]);
colorbar
colormap jet
title("Intensity at z = 0")
subplot(2,1,2);
plot(x,GB(:,ceil(L/2)));
title("Slice of Intensity at z = 0 & y = 0");
xlabel("x");
ylabel("I/Io");
ylim([0 1]);
xlim([-5*wo 5*wo]);
%------------------------------------------
%Part 2
%------------------------------------------
%n = 2^nextpow2(L);
Fs = 1/step;
k_step = Fs/L;
kx = linspace(-Fs/2,Fs/2,L)*2*pi;
ky = linspace(-Fs/2,Fs/2,L)*2*pi;
%kx = 2*pi*[-Fs/2:k_step:0+k_step k_step:k_step:Fs/2+k_step];
%ky = 2*pi*[-Fs/2:k_step:0+k_step k_step:k_step:Fs/2+k_step];
%kx = 2*pi*(1./x);
%ky = 2*pi*(1./y);
[kxx,kyy] = meshgrid(kx,ky);
Ukxky = fftshift(fft2(U(X,Y)));                                                                                         
%[FX,FY] = meshgrid(Fv,Fv);
%imagesc(kx,ky,abs(Ukxky)/n)
%colorbar
%caxis([0 1]);
%subplot(2,1,2);
%imagesc(kx,ky,angle(Ukxky))
%colorbar
kz = k - ((kxx.^2+kyy.^2)./(2*k));
%kz = (k^2 - kxx.^2-kyy.^2).*0.5;
U_z = @(z) abs(ifft2(ifftshift(Ukxky.*exp(-1i.*z.*kz)))).^2;
figure(2)
z_vals = [0.5 1 2];
for plotId = 1 : 3
    subplot(3,2,2*plotId-1);
    z_val = z_vals(plotId)*zo;
    Uz = U_z(z_val);
    imagesc(x,y,Uz)
    format short
    title(sprintf("Intensity at z = %.1f zo",z_vals(plotId)));
    xlabel("x");
    ylabel("y");
    colorbar
    caxis([0 1])
    colormap jet
    subplot(3,2,2*plotId);
    plot(x,Uz(:,ceil(L/2)));
    title(sprintf("Slice of Intensity at z = %.1f zo & y = 0",z_vals(plotId)));
    xlabel("x");
    ylabel("I/Io");
    ylim([0 1]);
    xlim([-5*wo 5*wo]);
end