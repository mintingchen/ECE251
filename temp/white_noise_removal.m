l = 10; % second
% how the white noise signal is generated?
% how to calculate the SNR?
% Perform 1-10 iteration and plot the output
% may need to record the clean sound again
file = "input/breath_clean.m4a";
samples = [l*48000,2*l*48000-1];
[x_data,fs] = audioread(file,samples);
t = 0:1/fs:l-1/fs;
x_data= x_data(:,1)*5;
%sound(x_data,fs)

lev = 10; 
snr_list = zeros(10,1);

for i = 1:lev
    y = awgn(x_data,5,'measured'); % Add white Gaussian noise. figure(1)
    %plot(t,y,t,x_data)
    wname = 'coif5'; 
    tree = wpdec(y,lev,wname); 
    det1 = wpcoef(tree,2014);
    sigma = median(abs(det1))/0.6745; 
    alpha = 1.8;
    thr = wpbmpen(tree,sigma,alpha) 
    keepapp = 1;
    xd = wpdencmp(tree,'s','nobest',thr,keepapp);
    noise = xd-x_data;
    snr_list(i)= 10*log10(sum(x_data.^2)/sum(noise.^2));
end
%%
noise_xy = y-x_data;
ori_snr = 10*log10(sum(x_data.^2)/sum(noise_xy.^2))
%%
index = linspace(1,10,10)
plot(index, snr_list)

% figure
% plot(t,x_data,t,xd,'g') 
% xlabel('time') 
% ylabel('Amplitude') 
% title('signal denoise');



