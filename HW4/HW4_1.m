% Shabab Ahmed 

clear all; close all; clc 

% loading the images into Matlab 

cropped = dir('CroppedYale'); % directory 
cropped_f = cropped(4:end); 

images = [];
for i = 1:length(cropped_f)
    if i<9
        strName = ['yaleB0' num2str(i)];
    else
        strName = ['yaleB' num2str(i)];
    end
    sub = dir(['CroppedYale', '/' ,strName]); 
    sub_file = sub(4:end); 
    for i = 1:length(sub_file)
        file = (sub_file(i).name); % get attribute name from struct 
        image1 = imread(file);  % read the file 
        sz= size(image1);
        sz1 = sz(1);
        sz2 = sz(2);
        sz = sz1*sz2;
        image = reshape(image1, sz,1); % reshaping image 
        images = [images, image];
    end
end



images =  double(images); 
[U,S,V] =svd(images,'econ'); 
sig = diag(S); 



figure(1)

subplot(2,2,1)
plot(sig, 'ko','Linewidth',[1.5])
title('Spectrum of Singular Values')
set(gca, 'Fontsize', [13])

subplot(2,2,2)
semilogy(sig,'ko','Linewidth',[1.5])
ylim([10^(2); 10^(6)]);
set(gca,'Fontsize',[13],'Ytick', [10^(2) 10^(3) 10^(4),...
10^(5), 10^(6)])
title('Semilogy plot of Singular Values')

subplot(2,1,2)
plot((sig/sum(sig))*100, '-.or')
xlabel('Number of modes')
ylabel('Percentage of energy captured')
ylim([0, 15]); 
title('Percentage of Energy at each mode')

figure(2)
title('First four modes') 
set(gca, 'Fontsize', [16])
% eigenfaces
for j = 1:4
    subplot(2,2,j)
    ef = reshape(U(:,j),sz1, sz2);
    pcolor(ef), axis off, shading interp, colormap(hot)
end
% 
m = randi([0 2269], 1, 3); 

figure(3)
title('Reconstruction of Faces')
for j = 1:length(m)
    i = m(j); 
    im = images(:,i); 
    im = reshape(im, sz1, sz2); 
    im = uint8(im); 
    subplot(3,2,(2*j)-1)
    imshow(im);
    A = U*S(:,1:50)*V(:,1:50)'; 
    A = A(:,i); 
    A = reshape(A, sz1, sz2); 
    A = uint8(A);
    subplot(3,2,2*j)
    imshow(A)
    
end

%------------------------------------------------------------------------

uncropped = dir('yalefaces_uncropped'); % directory
uncropped_f = uncropped(end);  

images_u = [];
sub = [uncropped_f.name]; % get attribute name from struct
subfile =dir(['yalefaces_uncropped', '/', sub]);
sub_file = subfile(4:end); 
for i = 1:length(sub_file)
        file = (sub_file(i).name); % get attribute name from struct 
        image1 = imread(file);  % read the file 
        sz= size(image1);
        sz1 = sz(1);
        sz2 = sz(2);
        sz = sz1*sz2;
        image = reshape(image1, sz,1); % reshaping image 
        images_u = [images_u, image];
end
 

images_u =  double(images_u); 
[U,S,V] =svd(images_u,'econ'); 
sig = diag(S); 


figure(4)

subplot(2,2,1)
plot(sig, 'ko','Linewidth',[1.5])
title('Spectrum of Singular Values')
set(gca, 'Fontsize', [13])

subplot(2,2,2)
semilogy(sig,'ko','Linewidth',[1.5])
% ylim([10^(2); 10^(6)]);
% set(gca,'Fontsize',[13],'Ytick', [10^(2) 10^(3) 10^(4),...
% 10^(5), 10^(6)])
title('Semilogy plot of Singular Values')

subplot(2,1,2)
plot((sig/sum(sig))*100, '-.or')
xlabel('Number of modes')
ylabel('Percentage of energy captured')
% ylim([0, 15]); 
title('Percentage of Energy at each mode')

figure(5)
title('First four modes') 
set(gca, 'Fontsize', [16])


% eigenfaces 
modes = [2,4,6,8];
for i = 1:length(modes)
    subplot(2,2,i)
    k = modes(i)
    ef = reshape(U(:,k),sz1, sz2);
    pcolor(ef), axis off, shading interp, colormap(hot)
end
% 
m = randi([0 165], 1, 3); 

figure(6)
title('Reconstruction of Faces')
for j = 1:length(m)
    i = m(j); 
    im = images_u(:,i); 
    im = reshape(im, sz1, sz2); 
    im = uint8(im); 
    subplot(3,2,(2*j)-1)
    imshow(im);
    A_u = U*S(:,1:4)*V(:,1:4)'; 
    A_u = A_u(:,i); 
    A_u = reshape(A_u, sz1, sz2); 
    A_u = uint8(A_u);
    subplot(3,2,2*j)
    imshow(A_u)
    
end



