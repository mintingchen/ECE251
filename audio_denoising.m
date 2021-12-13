clear all
% read input files
l = 10; % second
samples = [l*48000*1.5,2.5*l*48000-1];

clean_file = 'input/clean.m4a';
[clean,fs] = audioread(clean_file,samples);
clean= clean(:,1)';
t = 0:1/fs:l-1/fs;

% input1: music
music_file = 'input/music.m4a';
[music_noise,fs] = audioread(music_file,samples);
music_noise = music_noise(:,1)';
scaler = sqrt(sum(clean.^2)/(10^0.5)/sum(music_noise.^2));
music_noise1 = music_noise*scaler;
music_input = clean + music_noise1; % Mixed clean breath with noise
music_snr_before = snr(clean,music_noise1)
%sound(movie_input,fs)

% input2: movie
movie_file = 'input/movie.m4a';
[movie_noise,fs] = audioread(movie_file,samples);
movie_noise = movie_noise(:,1)';
scaler = sqrt(sum(clean.^2)/(10^0.5)/sum(movie_noise.^2));
movie_noise1 = movie_noise*scaler;
movie_input = clean + movie_noise1; % Mixed clean breath with noise
movie_snr_before = snr(clean,movie_noise1)
%sound(movie_input,fs)

% input3: synthesis white gaussian noise
wnoise_snr_before = 5
wnoise_input = awgn(clean,wnoise_snr_before,'measured');
%sound(wnoise_input, fs);
%%
% Experiment 1: iteration 
wname = 'haar';
iter = 10;
music_snr_iter = iteration_snr(music_input, clean, wname, iter, music_snr_before,"music");
movie_snr_iter = iteration_snr(movie_input, clean, wname, iter, movie_snr_before,"movie");
wnoise_snr_iter = iteration_snr(wnoise_input, clean, wname, iter, wnoise_snr_before,"white");

figure
plot(0:iter, music_snr_iter,0:iter, movie_snr_iter,0:iter, wnoise_snr_iter);
% text(0:iter,music_snr_iter,num2str(music_snr_iter));
% text(0:iter,movie_snr_iter,num2str(movie_snr_iter));
title('Number of iterations vs SNR of denoised signal')
xlabel('Number of iterations') 
ylabel('SNR') 
legend('music','movie','white noise')
ylim([0 15])
%%
% Experiment 2: different type of wavelet
wavelet_list = ["coif5";"db5";"fk8";"haar";"sym6"];
music_snr_wavelet = wavelet_snr(music_input, clean, wavelet_list, 7, music_snr_before,"music");
movie_snr_wavelet = wavelet_snr(movie_input, clean, wavelet_list, 9, movie_snr_before,"movie");
wnoise_snr_wavelet = wavelet_snr(wnoise_input, clean, wavelet_list, 4, wnoise_snr_before,"white");
figure
plot(0:5, music_snr_wavelet,0:5, movie_snr_wavelet ,0:5, wnoise_snr_wavelet);
% text(0:5,music_snr_iter,num2str(music_snr_iter));
% text(0:iter,movie_snr_iter,num2str(movie_snr_iter));
title('Wavelet Choises vs SNR of denoised signal')
xlabel('Wavelet type') 
ylabel('SNR') 
legend('music','movie','white noise')
ylim([0 18])
xticks([0 1 2 3 4 5])
xticklabels({"original","coif5","db5","fk8","haar","sym6"})
%%
% Experiment 3: different initial SNR with best iteration and wavelet type
initial_snr = [-15;-10; -5; 0; 5; 10;15;20]
music_initial_snr = diff_initial_snr(music_noise, clean, "haar", 7, initial_snr,music_snr_before,"music");
movie_initial_snr = diff_initial_snr(movie_noise, clean, "haar", 9, initial_snr,movie_snr_before,"movie");
wnoise_initial_snr = diff_initial_snr(movie_noise, clean, "coif5", 4, initial_snr,wnoise_snr_before,"white");

a = size(initial_snr,1)-1;
figure
plot(0:a, music_initial_snr,0:a, movie_initial_snr ,0:a, wnoise_initial_snr,0:a,initial_snr','--');
% text(0:5,music_snr_iter,num2str(music_snr_iter));
% text(0:iter,movie_snr_iter,num2str(movie_snr_iter));
title('Initial SNR vs SNR of denoised signal')
xlabel('Initial SNR') 
ylabel('Output SNR') 
legend('music','movie','white noise')
ylim([-15 25])
xticks([0 1 2 3 4 5 6 7])
xticklabels({"-15","-10","-5","0","5","10","15","20"})

function xd = wavelet_denoising(input,lev,wname,noise_type)
tree = wpdec(input,lev,wname); 
% figure 
% plot(tree)
if noise_type == "music"
    index = 5;
    alpha = 1;
elseif noise_type == "movie"
    index = 5;
    alpha = 1;
else
    index = 1;
    alpha = 2;
end
det1 = wpcoef(tree,2^(lev));
sigma = median(abs(det1))/0.6745/index; %music/5 movie/10
thr = wpbmpen(tree,sigma,alpha);
keepapp = 1;
xd = wpdencmp(tree,'s','nobest',thr,keepapp);
end

function snr_iter = iteration_snr(input, clean, wname, iter, snr_before,noise_type)
    snr_iter = zeros(iter,1);
    % denoising algorithm
    for i = 1:iter
        xd = wavelet_denoising(input,i,wname,noise_type);
        noise_after = abs(xd - clean);
        snr_after = 10*log10(sum(clean.^2)/sum(noise_after.^2));
        snr_iter(i) = snr_after;
        if noise_type == "music"
            index = 7;
        elseif noise_type == "movie"
            index = 9;
        else
            index = 4;
        end

        if i == index
            fs = 48000;
            t = 0:1/fs:10-1/fs;
            %sound(xd*2,48000);
            figure
            subplot(2,1,1)
            plot(t, input, t, clean);
            legend('input','clean')
            xlabel('time') 
            ylabel('Magnitude') 
            title(noise_type+': clean signal vs input signal')
            subplot(2,1,2)
            plot(t, input, t, xd);
            legend('input','output')
            xlabel('time') 
            ylabel('Magnitude') 
            title(noise_type+': denoised signal vs input signal')
            %sound(xd,fs)
        end

    end
    snr_iter = [snr_before; snr_iter];
end

function snr_iter = wavelet_snr(input, clean, wavelet_list, lev, snr_before,noise_type)
    len_wavlist = size(wavelet_list,1)
    snr_iter = zeros(len_wavlist,1);
    % denoising algorithm
    
    for i = 1:len_wavlist
        xd = wavelet_denoising(input,lev,wavelet_list(i),noise_type);
        noise_after = abs(xd - clean);
        snr_after = 10*log10(sum(clean.^2)/sum(noise_after.^2));
        snr_iter(i) = snr_after;
        if noise_type == "music"
            index = 7;
        elseif noise_type == "movie"
            index = 9;
        else
            index = 4;
        end
        if i == index
            fs = 48000;
            t = 0:1/fs:10-1/fs;
            %sound(xd*2,48000);
            figure
            subplot(2,1,1)
            plot(t, input, t, clean);
            legend('input','clean')
            xlabel('time') 
            ylabel('Magnitude') 
            title('clean signal vs input signal')
            subplot(2,1,2)
            plot(t, input, t, xd);
            legend('input','output')
            xlabel('time') 
            ylabel('Magnitude') 
            title('denoised signal vs input signal')
            sound(xd,fs)
        end

    end
    snr_iter = [snr_before; snr_iter];
    
end

function snr_iter = diff_initial_snr(noise, clean, wavelet, lev, initial_snr_list,snr_before,noise_type)
    len_snrlist = size(initial_snr_list,1)
    snr_iter = zeros(len_snrlist,1);
    % denoising algorithm
    for i = 1:len_snrlist
        if noise_type == "white"
            input = awgn(clean,initial_snr_list(i),'measured');
        else
            scaler = sqrt(sum(clean.^2)/(10^(initial_snr_list(i)/10))/sum(noise.^2));
            input = clean+noise*scaler;
            snr(clean,noise*scaler)
        end
        xd = wavelet_denoising(input,lev,wavelet,noise_type);
        noise_after = abs(xd - clean);
        snr_after = 10*log10(sum(clean.^2)/sum(noise_after.^2));
        snr_iter(i) = snr_after;
        if noise_type == "music"
            index = 7;
        elseif noise_type == "movie"
            index = 9;
        else
            index = 4;
        end
        if i == index
            fs = 48000;
            t = 0:1/fs:10-1/fs;
            %sound(xd*2,48000);
            figure
            subplot(2,1,1)
            plot(t, input, t, clean);
            legend('input','clean')
            xlabel('time') 
            ylabel('Magnitude') 
            title('clean signal vs input signal')
            subplot(2,1,2)
            plot(t, input, t, xd);
            legend('input','output')
            xlabel('time') 
            ylabel('Magnitude') 
            title('denoised signal vs input signal')
            sound(xd,fs)
        end

    end
%     snr_iter = [snr_before; snr_iter];  
end

