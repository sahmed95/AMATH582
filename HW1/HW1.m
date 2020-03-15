clear all; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes 2^6
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1];%making it 2pi periodic
ks=fftshift(k); %shifting


%creating mesh grids
[X,Y,Z]=meshgrid(x,y,z); 
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

%Part 1: Averaging to find the frequency signature

%initializing for averaging
uave = zeros;
for j =1:20
    uave = uave + fftn(reshape(Undata(j,:),n,n,n));  %adding the signal
end
utaves =fftshift(uave)/20;  %averaging the signal

figure(1)
close all, isosurface(Kx,Ky,Kz,abs(utaves)/max(abs(utaves(:))),0.6)
axis([-6 6 -6 6 -6 6]), grid on, drawnow
set(gca, 'FontSize', 20);
xlabel("Kx"); ylabel("Ky"); zlabel("Kz");


%finding the maximum and the corresponding index:
[m,index] = max(abs(utaves(:)));
[i_x,i_y,i_z] = ind2sub(size(abs(utaves)), index);

% finding the corresponding wavenumbers (frequency signature)
k_x = Kx(i_x,i_y,i_z); 
k_y = Ky(i_x, i_y, i_z); 
k_z = Kz(i_x,i_y,i_z);

freq_signature = [k_x, k_y, k_z]

%2. Filtering the data to denoise and find trajectory
tau =0.2;

% Using the center frequency to create the filter
filter = exp(-tau*((Kx - k_x).^2 + (Ky - k_y).^2 + (Kz - k_z).^2));

% Initializing vectors for the position
X_pos = [];
Y_pos = [];
Z_pos = [];
for j = 1:20
    un(:,:,:) = reshape(Undata(j,:),n,n,n);
    Ut = fftn(un);  % Fourier transform
    Uts = filter.*fftshift(Ut); % filtering
    Un = ifftn(ifftshift(Uts));  % Inverse transform
    [m,index] = max(abs(Un(:))); % Finding the max
    [x_p, y_p, z_p] = ind2sub(size(abs(Un)),index);
    X_pos = [X_pos; X(x_p,y_p, z_p)];
    Y_pos = [Y_pos; Y(x_p, y_p, z_p)];
    Z_pos = [Z_pos; Z(x_p, y_p, z_p)];
end

%Plotting isosurface of inverse Fourier transform of last measurement
figure(2)
close all, isosurface(X,Y,Z, abs(Un)/max(abs(Un(:))),0.6)
axis([-L L -L L -L L]), grid on, drawnow
set(gca, 'FontSize', 20);
xlabel("X"); ylabel("Y"); zlabel("Z");

%Plotting the trajectory
figure(3)
plot3(X_pos, Y_pos, Z_pos)
grid('on')
set(gca, 'FontSize', 20);
xlabel("X"); ylabel("Y"); zlabel("Z");


%Acoustic wave location(final position of marble): 
position = [X_pos(20), Y_pos(20), Z_pos(20)];



