% Shabab Ahmed
% HW 2 582
% Part 1

%--------------------------------------------------------------------
clear all; close all; clc

% load signal

load handel

% take transpose of signal and divide by 2

v = y'/2;



% make signal even length
v= v(1:end-1);
n = length(v);
L = (n)/Fs; %final time 
t = linspace(0,L,n); %time domain
k = (2*pi/L)*[0:n/2-1 -n/2:-1];  % make k 2pi periodic
ks = fftshift(k);                % shift 

%--------------------------------------------------------------------

% Plotting signal of interest 

figure(1)
subplot(2,1,1) % Time domain 
plot((1:length(v))/Fs,v, 'k')
xlim([0 t(end-1)])
set(gca, 'Fontsize', [20]),
xlabel('Time (t)'), ylabel('v(t)')
title('Signal of interest (v)')

% Fourier transform of signal 
vt = fft(v); 

subplot(2,1,2) %Fourier domain 
plot(ks, abs(fftshift(vt))/max(abs(vt)),'k');
set(gca, 'Fontsize', [20])
xlabel('Frequency (k)');
ylabel('FFT(v)');
title('FFT of v')

%--------------------------------------------------------------------

% Gabor filtering using Gaussian window 

figure(2)


a = 500;
vgt_spec=[];
dt = 0.1;
tslide=0:dt:t(end-1);
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); % Gabor window
    vg=g.*v;                    % product of signal and filter
    vgt=fft(vg);                
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
%     subplot(3,1,1), plot(t,v,'k',t,g,'r')
%     xlabel('Time(sec)'); ylabel('Amplitude')
%     title('Signal with Gaussian Filter')
% 
%     subplot(3,1,2), plot(t,vg,'k')
%     xlabel('Time(sec)'); ylabel('Amplitdue')
%     title('Product of Signal and Filter')
%     subplot(3,1,3), plot(ks,abs(fftshift(vgt)/max(abs(fftshift(vgt)))))
%     xlabel('Time(sec)'); ylabel('Amplitude')
%     title('FFT of Product')
%     axis([-50 50 0 0.3])
%     drawnow
%     pause(0.1)
end
figure(3)
subplot(2,2,1)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)'); ylabel('Frequency (\omega)')
title('a = 500')
colormap(hot)


%--------------------------------------------------------------------

vgt_spec=[];
a = 250;
tslide=0:dt:t(end-1);
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); % Gabor
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    
end

subplot(2,2,2)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)'); ylabel('Frequency (\omega)')
title('a = 250')
colormap(hot)

%--------------------------------------------------------------------

vgt_spec=[];
a = 50; 
tslide=0:dt:t(end-1);
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); % Gabor
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    
end

subplot(2,2,3)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)'); ylabel('Frequency (\omega)')
title('a = 50')
colormap(hot)

%--------------------------------------------------------------------


vgt_spec=[];
a = 2; 
tslide=0:dt:t(end-1);
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); % Gabor filter Gaussian
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    
end

subplot(2,2,4)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)'); ylabel('Frequency (\omega)')
title('a = 2')
colormap(hot)

%--------------------------------------------------------------------


% Undersampling

figure(4)
a = 250; 
dt = 0.5;
vgt_spec=[];
tslide=0:dt:t(end-1);
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); % Gabor filter Gaussian
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
   
end

subplot(2,1,1)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)'); ylabel('Frequency (\omega)')
title('Undersampled Signal, a = 250, dt = 0.5')
colormap(hot)

% Oversampling 
a = 50; 
dt = 0.01;
vgt_spec=[];
tslide=0:dt:t(end-1);
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); % Gabor filter
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
   
end

subplot(2,1,2)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)'); ylabel('Frequency (\omega)')
title('Oversampled Signal, a = 250, dt = 0.01')
colormap(hot)
%--------------------------------------------------------------------


% Mexican hat wavelet 


vgt_spec=[];
a = 0.05; 

% function for mexican hat wavelet 

m = @(x) (1-(x/a).^2).*exp(-(x/((sqrt(2)*a))).^2).*(1/sqrt(a));
dt = 0.1;
tslide=0:dt:t(end-1);
for j=1:length(tslide)
    x = t-tslide(j);
    g=m(x); % Mexican hat
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];

end

figure(5)
y = t-tslide(50); 
g = m(y); 
subplot(2,1,1)
plot(t,v,'k',t,g,'r')
xlabel('Time(sec)'); ylabel('Amplitude')
axis([4 6 -3 6])
title('Signal with Mexican hat wavelet')

subplot(2,1,2)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)'); ylabel('Frequency (\omega)')
title('Mexican hat wavelet Spectogram, a =0.05, dt = 0.1')
colormap(hot)

%-------------------------------------------------------------------------

% Shannon window 
a = 3;

%function for Shannon window

s = @(x) (a*abs(x) <= 0.5);

vgt_spec = [];
dt = 0.1;
tslide=0:dt:t(end-1);
for j=1:length(tslide)
    y = t-tslide(j);
    g=s(y);  
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    
end

figure(6)
y = t-tslide(40); 
g = s(y); 
subplot(2,1,1)
plot(t,v,'k',t,g,'r')
xlabel('Time(sec)'); ylabel('Amplitude')
set(gca, 'Fontsize', [15]),
xlim([0, t(end)])
title('Signal with Shannon filter')

subplot(2,1,2)
pcolor(tslide,ks,vgt_spec.'), shading interp

xlabel('Time(sec)'); ylabel('Frequency (\omega)')
title('Shanon Window Spectogram, a =3, dt = 0.1')
set(gca, 'Fontsize', [15]),
colormap(hot)