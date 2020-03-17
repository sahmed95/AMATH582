%Part 2
%Test 3

clear all; close all; clc
%%
x = audioread('Alternative_rock.wav');
y = audioread('Grunge.wav');
z = audioread('Classical.wav'); 
tr = 75; 

% audio normalization 
xsum = sum(x,2);
mLx = max(abs(x(:,1)));
mRx = max(abs(x(:,2)));
max_x = max([mLx mRx]);
x = xsum*max_x;

ysum = sum(y,2);
mLy = max(abs(y(:,1)));
mRy = max(abs(y(:,2)));
max_y = max([mLy mRy]);
y = ysum*max_y;

zsum = sum(z,2);
peakAz = max(abs(zsum));
peakLz = max(abs(z(:,1)));
peakRz = max(abs(z(:,2)));
max_z = max([peakLz peakRz]);
z = ysum*max_z;

% 
Fs = length(x); 


x= x'/2; 
y= y'/2;
z = z'/2;

r = mod(length(x),2);

% make even length
if  r == 1
    x = x(1:end-1); 
    y = y(1:end-1);
    z = z(1:end-1);
    
end


n = length(x); 
L = tr;    % final time 
% discretized time vectors

t = linspace(0,L, n); 

k = (2*pi/L)*[0:n/2-1 -n/2:-1];  %rescale 
ks = fftshift(k); %shift
x=single(x); y=single(y); z =single(z);
%%
a = 250; % window width
dt = 1.5; 
xgt_spec = []; ygt_spec = []; zgt_spec = [];
xg_spec =[]; yg_spec= []; zg_spec =[];
tslide = 0:dt:t(end-1);

for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); % Gabor
    xg = g.*x; xgt =fft(xg); 
    yg=g.*y; ygt=fft(yg);
    zg = g.*z; zgt = fft(zg);
    xgt =single(xgt); ygt =single(ygt); zgt=single(zgt);
    xgt_spec=[xgt_spec; abs(fftshift(xgt))];
    ygt_spec=[ygt_spec; abs(fftshift(ygt))];
    zgt_spec=[zgt_spec; abs(fftshift(zgt))];
    
    % make sure each row in spectogram is a song
    % if time is a multiple of 5s we have a new song
    if mod(j+5,5)== 0
        xg_spec = [xg_spec; xgt_spec];
        yg_spec = [yg_spec; xgt_spec];
        zg_spec = [zg_spec; xgt_spec];
        xgt_spec = [];
        ygt_spec = [];
        zgt_spec = []; 
    end 
end

Spect = [xg_spec; yg_spec; zg_spec];
%%
[U,S,V] = svd(Spect, 'econ');
sig = diag(S); 
plot(sig/sum(sig)*100,  '-.or')

figure(2)
plot(sig, 'ko')
title('Spectrum of Singular Values')
set(gca, 'Fontsize', [14])



%%
features = 5;

train_s = 10;
test_s = 5; 

% 
figure(3)
answer = [ones((test_s),1); 2*ones(test_s,1); 3*ones(test_s,1)]; 
bar(answer)
title('Correct Classification')
xlim([0.5 15.5])

% initializing vectors for accuracy
acc_LDA = [];
acc_nb = [];

trials =500; %number of trials
for j = 1:trials
    % randomly choosing songs
    q = linspace(1,15,15);
    q1 = randperm(15); 
    q2 = randperm(15);
    q3 = randperm(15) ;
    q1 = q1(1:train_s);
    q11 = setdiff(q,q1);
    q2 = q2(1:train_s)+15; 
    q=linspace(16,30,15);
    q22 = setdiff(q,q2);
    q3 = q3(1:train_s)+30;
    q = linspace(31,45,15);
    q33 = setdiff(q,q3);
    
    % create training set 
    
    xtrain = [];ytrain =[]; ztrain=[];
    for i = 1:length(q1)
        sl = 80000;
        xt = V(q1(i), 1:features);
        xtrain = [xtrain; xt]; 
        yt = V(q2(i), 1:features);
        ytrain = [ytrain; yt];
        zt = V(q3(i), 1:features);
        ztrain = [ztrain; zt];
    end 
    % create test set 
    
    xtest = []; ytest =[]; ztest = [];
    for k =1:length(q11)
        xt = V(q11(k), 1:features);
        xtest = [xtest; xt]; 
        yt = V(q22(k), 1:features);
        ytest = [ytest; yt];
        zt = V(q22(k), 1:features);
        ztest = [ztest; zt]; 
    end 
    test = [xtest; ytest; ztest];
    train = [xtrain; ytrain; ztrain];
    
    % grouping 
    ctrain = [ones(train_s,1); 2*ones(train_s,1); 3*ones(train_s,1)];
    
    % Linear Discriminant Analysis 
    
    pre =classify(test, train, ctrain);
    % accuracy
    acc = 0;
    for l= 1:length(pre)
        if pre(l) == answer(l)
            acc =acc+1; 
        end
    end 
    acc_LDA = [acc_LDA;acc/length(pre)*100];
    
    % Naive Bayes 
    nb = fitcnb(train, ctrain);
    pre1 = nb.predict(test); 
    
    % accuracy 
    
    acc_1 = 0;
    for p = 1:length(pre1)
        if pre1(p) == answer(p)
            acc_1=acc_1 +1;
        end
    end
    acc_nb = [acc_nb;acc_1/length(pre1)*100];
    
end

% Average accuracy 

accuracy_LDA = sum(acc_LDA)/length(acc_LDA); 
accuracy_nb = sum(acc_nb)/length(acc_nb);


% 
figure(4)
bar(pre);
title('Classification of Test set (LDA)')
set(gca,'Fontsize', [14])
ylim([0 3])
xlim([0.5 15.5])
% 
figure(5)
bar(pre1)
title('Classification of Test set (Naive Bayes)')
set(gca,'Fontsize', [14])
ylim([0 3])
xlim([0.5 15.5])

figure(6)
bar(ctrain)
title('Classification of Training set')
set(gca,'Fontsize', [14])
ylim([0 3])
xlim([0.5 30.5])

