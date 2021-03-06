% Noisy case 

clear all; close all; clc 

for k = 1:3
    load(['cam' num2str(k) '_2'])
end
vidFrames1 = vidFrames1_2;
vidFrames2 = vidFrames2_2;
vidFrames3 = vidFrames3_2;

%[height, width, [R,G,B], number of frames]

[h1, w1, c1,n1] = size(vidFrames1);
[h2, w2, c2,n2] = size(vidFrames2);
[h3, w3, c3,n3] = size(vidFrames3);


% obtain number of frames for the cameras
numFrames1 = n1;
numFrames2 = n2;
numFrames3 = n3;

numFrames = min([numFrames1, numFrames2, numFrames3]);

for k = 1 : numFrames1
    mov1(k).cdata = vidFrames1(:,:,:,k);
    mov1(k).colormap = [];
end

for k = 1 : numFrames2
    mov2(k).cdata = vidFrames2(:,:,:,k);
    mov2(k).colormap = [];
end

for k = 1 : numFrames3
    mov3(k).cdata = vidFrames3(:,:,:,k);
    mov3(k).colormap = [];
end


% Initialize position vectors for camera 1

X1 = [];
Y1 = [];

% window (box) size 
windowsize = 5;

for k=1:numFrames1
    gscale = frame2im(mov1(k));
    gscale = rgb2gray(gscale);  % convert to grayscale
    gscale(:,1:275) = 0;        % zero(black) out specific region in video
    maxavg = 0;                 % set maximum average pixel value of window
    pos = [];                    % initial position
    
    
    for j = 1:windowsize:h1
        for i = 1:windowsize:w1
            % if height and width out of bounds set appropriate window
            if (j+windowsize>h1) & (i+windowsize>w1)
                box = gscale(j:h1,i:w1);
            
            % if height out of bounds set appropriate window
            elseif (j+windowsize>h1)
                box = gscale(j:h1,i:i+windowsize);
            
            % if width out of bounds set appropriate window
            elseif (i + windowsize > w1)
                box = gscale(j:j+windowsize, i:w1);
                
            % normal window     
            else
                box = gscale(j:j+windowsize, i:i+windowsize);
            end
            % average the pixel value in the window 
            avg = sum(box(:));
            avg = avg/(windowsize*windowsize);            
            
            % update maximum average and position found for each window
            if avg > maxavg
                maxavg = avg;
                pos = [j,i];
             
           
            end
        end
    end
    
    % extract bottom left position of box with maximum average pixel value
    j= pos(1);
    i = pos(2);
    
    %recreate the window of maximum average pixel value
    subframe = gscale(j:j+windowsize,i:i+windowsize);
    
    [m, ind] = max(subframe(:)); % find maximum pixel intensity
    [x1 y1] = ind2sub(size(subframe), ind); % find corresponding indices
    X1 = [X1, x1+j];
    Y1 = [Y1, y1+i];
end

%----------------------------------------------------------------------
% Initialize position vectors for camera 2 

X2 = [];
Y2 = [];

% window (box) size

windowsize = 5;

for k=1:numFrames2
    gscale = frame2im(mov2(k));
    gscale = rgb2gray(gscale);  % convert to grayscale
    
    % Zero (black) out bright pixels that are not the paint can
    
    gscale(:,1:275) = 0;
    gscale(:,330:end) =0;
    gscale(1:200,:) = 0;
    
    maxavg = 0;
    pos = [];
    
    for j = 1:windowsize:h1
        for i = 1:windowsize:w1
            if (j+windowsize>h1) & (i+windowsize>w1)
                box = gscale(j:h1,i:w1);
            elseif (j+windowsize>h1)
                box = gscale(j:h1,i:i+windowsize);
            elseif (i + windowsize > w1)
                box = gscale(j:j+windowsize, i:w1);
            else
                box = gscale(j:j+windowsize, i:i+windowsize);
            end
              
            avg = sum(box(:));
            avg = avg/(windowsize*windowsize);            
            if avg > maxavg
                maxavg = avg;
                pos = [j,i];
             
           
            end
        end
    end
    j= pos(1);
    i = pos(2);
    subframe = gscale(j:j+windowsize,i:i+windowsize);
    [m, ind] = max(subframe(:)); % find maximum pixel intensity
    [x2 y2] = ind2sub(size(subframe), ind); % find corresponding indices
    X2 = [X2, x2+j];
    Y2 = [Y2, y2+i];
end

%--------------------------------------------------------------------
% Initialize position vectors for camera 3

X3 = [];
Y3 = [];

% window (box) size

windowsize = 5;

for k=1:numFrames3
    gscale = frame2im(mov3(k));
    gscale = rgb2gray(gscale);  % convert to grayscale
   
    % Zero (black) out bright pixels that are not the paint can
    
    gscale(1:150,:) = 0;
    gscale(:,1:230) = 0; 
    maxavg = 0;
    pos = [];
    for j = 1:windowsize:h1
        for i = 1:windowsize:w1
            if (j+windowsize>h1) & (i+windowsize>w1)
                box = gscale(j:h1,i:w1);
            elseif (j+windowsize>h1)
                box = gscale(j:h1,i:i+windowsize);
            elseif (i + windowsize > w1)
                box = gscale(j:j+windowsize, i:w1);
            else
                box = gscale(j:j+windowsize, i:i+windowsize);
            end
              
            avg = sum(box(:));
            avg = avg/(windowsize*windowsize);            
            if avg > maxavg
                maxavg = avg;
                pos = [j,i];
             
           
            end
        end
    end
    j= pos(1);
    i = pos(2);
    subframe = gscale(j:j+windowsize,i:i+windowsize);
    [m, ind] = max(subframe(:)); % find maximum pixel intensity
    [x3 y3] = ind2sub(size(subframe), ind); % find corresponding indices
    X3 = [X3, x3+j];
    Y3 = [Y3, y3+i];
end



% Aligning the position vectors so that the lowest point matches up 

% find the lowest position in the first 30 frames

[m1, i1] = min(X1(1:30)); 
[m2, i2] = min(X2(1:30)); 
[m3, i3] = min(X3(1:30));  


% align and crop the position vectors to 280 frames 

X1 = X1(i1:i1+280); Y1 = Y1(i1:i1+280);
X2 = X2(i2:i2+280); Y2 = Y2(i2:i2+280);
X3 = X3(i3:i3+280); Y3 = Y3(i3:i3+280); 


% Plotting the position vectors 

figure(2)
% Camera 1
subplot(3,2,1)
plot(Y1)
axis([0 210 0 600])
ylabel('Position in X')
xlabel('Frame number')
title('Camera 1')

subplot(3,2,2)
plot(X1)
axis([0 210 0 600])
ylabel('Position in Y')
xlabel('Frame number')
title('Camera 1')

% Camera 2
subplot(3,2,3)
plot(Y2) 
axis([0 210 0 600])
ylabel('Position in X')
xlabel('Frame number')
title('Camera 2')

subplot(3,2,4)
plot(X2) 
axis([0 210 0 600])
ylabel('Position in Y')
xlabel('Frame number')
title('Camera 2')

% Camera 3
subplot(3,2,5)
plot(Y3) 
axis([0 210 0 600])
ylabel('Position in X')
xlabel('Frame number')
title('Camera 3')

subplot(3,2,6)
plot(X3)
axis([0 210 0 600])
ylabel('Position in Y')
xlabel('Frame number')
title('Camera 3')



X = [X1;Y1;X2;Y2;X3;Y3]; % creating the data matrix 
[m,n] = size(X);  % compute data size
mn = mean(X,2); % mean for each row
X = X - repmat(mn,1,n); % subtract mean 

[u,s,v] = svd(X'/sqrt(n-1),0); % perform the svd
lambda = diag(s).^2;  % produce diagonal variances

sig =diag(s); % singular values 



% Plot percentage variance explained by each mode 
figure(3) 
plot(lambda/sum(lambda)*100, '-.or')
axis([1 6 0 100])
set(gca,'Fontsize',[13],'Ytick',[0 10 20 30 40 50 60 70 80 90 100],...
'Xtick',[1 2 3 4 5 6]);
ylabel('% of variance explained')
xlabel('Modes')
title('Percentage variance explained by each mode')



% Plot singular values on a standard and log plot
figure(4)

% Standard plot
subplot(3,2,1) 
plot(sig, 'ko','Linewidth',[1.5])
axis([0 10 0 100])
set(gca,'Fontsize',[10],'Xtick',[0 5 10])
text(8,50,'(a)','Fontsize',[13])
title('Singular values on a standard plot')

% Log plot
subplot(3,2,2), semilogy(sig,'ko','Linewidth',[1.5])
axis([0 10 10^(-2) 10^3])
set(gca,'Fontsize',[10],'Ytick',[10^(-2) 10^(-1) 10^(0) 10^(1) 10^(2)... 
10^3],'Xtick',[0 5 10]);
text(8,8,'(b)','Fontsize',[13])
title('Singular values on a log plot')

% Plot vertical components of first two PCA modes
subplot(3,2,3)
plot([v(1,1),v(3,1), v(5,1)], 'bo', 'Linewidth', [2]); hold on;
plot([v(1,2), v(3,2),v(5,2)],'ro', 'Linewidth',[2])
legend('mode 1', 'mode 2','Location','NorthWest')
set(gca,'Fontsize',[10], 'Xtick', [1 2 3])
text(2.8,0.5,'(c)','Fontsize',[13])
xlabel('Camera number')
title('Vertical components of PCA modes')

% Plot horizontal components of first two PCA modes
subplot(3,2,4)
plot([v(2,1),v(4,1), v(6,1)], 'bo','Linewidth',[2]); hold on;
plot([v(2,2), v(4,2),v(6,2)],'ro', 'Linewidth',[2])
set(gca,'Fontsize',[10], 'Xtick', [1 2 3])
legend('mode 1', 'mode 2','Location','NorthWest')
text(2.8,0.3,'(d)','Fontsize',[13])
title('Horizontal components of PCA modes')
xlabel('Camera number')

% Plot of time evolution behavior of PCA modes

subplot(3,1,3)
plot(u(:,1), 'b', 'Linewidth',[2]); hold on
plot(u(:,2), 'r', 'Linewidth',[2]);
xlim([0 205])
text(190,0.15,'(e)','Fontsize',[13])
legend('mode 1', 'mode 2','Location','SouthEast')
title('Time evolution behavior of the important PCA modes')
xlabel('Frame number')

