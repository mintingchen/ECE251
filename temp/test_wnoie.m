clear all
% read input files
l = 10; % second
samples = [l*48000*1.5,2.5*l*48000-1];

clean_file = 'input/clean.m4a';
[clean,fs] = audioread(clean_file,samples);
clean= clean(:,1);
t = 0:1/fs:l-1/fs;

% input3: synthesis white gaussian noise
wnoise_snr_before = 5
wnoise_input = awgn(clean,wnoise_snr_before,'measured');
%sound(wnoise_input, fs);

% Experiment 1: iteration 
wname = 'haar';
iter = 10;
wnoise_snr_iter = iteration_snr(wnoise_input, clean, wname, iter, wnoise_snr_before,"white");

figure
plot(0:iter, wnoise_snr_iter);
% text(0:iter,music_snr_iter,num2str(music_snr_iter));
% text(0:iter,movie_snr_iter,num2str(movie_snr_iter));
xlabel('iter') 
ylabel('SNR') 
function xd = wavelet_denoising(input,lev,wname,noise_type)
tree = wpdec(input,lev,wname); 
% figure 
% plot(tree)
det1 = wpcoef(tree,2^(lev));
sigma = median(abs(det1))/0.6745; %music/5 movie/10
alpha =2;
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
            index = 7;
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

function snr_iter = wavelet_snr(input, clean, wavelet_list, lev, snr_before)
    len_wavlist = size(wavelet_list,1)
    snr_iter = zeros(len_wavlist,1);
    % denoising algorithm
    
    for i = 1:len_wavlist
        xd = wavelet_denoising(input,lev,wavelet_list(i));
        noise_after = abs(xd - clean);
        snr_after = 10*log10(sum(clean.^2)/sum(noise_after.^2));
        snr_iter(i) = snr_after;
        
        if i == 8
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
    figure
    plot(0:len_wavlist, snr_iter);
    te = text(0:len_wavlist,snr_iter,num2str(snr_iter));
    xlabel('wavelet') 
    ylabel('SNR') 
    title('SNR vs wavelet')
    xticks([0 1 2 3 4 5])
    xticklabels({"original","coif5","db5","fk8","haar","sym6"})

end
