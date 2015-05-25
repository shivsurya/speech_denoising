%Multi-band Spectral subtraction
%VAD used is a simple spectral threshold as described in paper.
%-Shiv Surya
%

close all
clear all

tic;
%read noisy data
[data,Fs]=audioread('noisy.wav');

%premphasis
data=filter([1 -0.96],1,data);

%noise silence period
len_ana=0.025; %window length of analysis is 25ms
init_sil=0.25;% initial silence in seconds
overlap_per=0.6; %percentage overlap in windows for enframing 

%Use hanning window as described in paper
win=hamming(len_ana*Fs);
win_len=len_ana*Fs;

%create a function that generates overlapping windows of the input signal
%taking into account the initial silence period
%what do I need?signal, overlap, windowlength,make this a function!!
%need the vector column-wise as fft is computed column-wise by matlab

win_overlap=floor(overlap_per*length(win));

%neglect the last few files less than win
noFr=floor((length(data)-length(win))/win_overlap+1);

%generate indices of the overlapping matrices
matInd=(repmat(1:length(win),noFr,1)+repmat((0:win_overlap:(noFr-1)*win_overlap)',1,length(win)))';
enframe_sigMat=data(matInd);

clear matInd;

%multiply the frame by the window
enframe_sigMat=enframe_sigMat.*repmat(win,1,noFr);
%clear win;


 FFT,compute mag,
% bias, halfwave, reduce residual and them VAD using spectral subtraction
% need to keep updating noise spectral average estimate as I detect using
% VAD

%compute FFT using  the same FFT length as window length
sigfftMat=fft(enframe_sigMat,win_len);
%fft is symmetric so compute magnitude of half the length of each column+DC

sigfftmag = abs(sigfftMat(1:floor(size(sigfftMat,1)/2)+1,:)).^2;
sigfftphase= angle(sigfftMat(1:floor(size(sigfftMat,1)/2)+1,:));


%compute the noise average from the initial silence segments
no_iniSilFr=floor(init_sil*Fs/win_overlap-1);

%compute mean rowwise for the intial silence period. This can be updated
%when the VAD detects the period. To keep track of the length of the
%silence period to update it, I use sil_len

noiseAvg=mean(sigfftmag(:,1:no_iniSilFr),2);
noise_len=no_iniSilFr;




%Intialize variables  for denoising(notation follows the paper )
upp_snr=20;
low_snr=-5;
upp_alpha=3;
low_alpha=6;
s=(low_alpha-upp_alpha)/(low_snr-upp_snr);
alpha0=upp_alpha-s*upp_snr;

%initialize delta(fill bottom down to avoid any missing value due to round off) 
del=1.5*ones(floor(win_len/2)+1,1);
del(1:floor((Fs/2.0-2000)*win_len/Fs))=2.5;
del(1:floor(1000*win_len/Fs))=1;
band_ind=[floor(1000*win_len/Fs) floor((Fs/2.0-2000)*win_len/Fs) ];

Gain=0.002;
dbthresh=6;
temp_thresh=8;
counter=9;
isspeech=0;

sig_recon=zeros(floor(win_len/2)+1,noFr);

%process complete signal frame wise

for i=1:noFr
    %VAD
    [isspeech,counter]=spectralDistVad(sigfftmag(:,i).^0.5,noiseAvg.^0.5,counter,dbthresh,temp_thresh);
    if(isspeech==0)
    %update estimate of average noise h
        noise_len=noise_len+1;
        noiseAvg=((noise_len-1)*noiseAvg+sigfftmag(:,i))./noise_len;
    end
    %compute snr and calculate dynamic value of alpha
%     snr1=10*log(sigfftmag(1:band_ind(1),i)./noiseAvg(1:band_ind(1)));
%     snr2=10*log(sigfftmag(band_ind(1)+1:band_ind(2),i)./noiseAvg(band_ind(1)+1:band_ind(2)));
%     snr3=10*log(sigfftmag(band_ind(2)+1:end,i)./noiseAvg(band_ind(2)+1:end));
%     snr=[snr1;snr2;snr3];
    snr=10*log(sigfftmag(:,i)./noiseAvg);
    alpha=alpha0-snr*s;
    alpha(snr>upp_snr)=upp_alpha;
    alpha(snr<low_snr)=low_alpha;
    sig_recon(:,i)=max((sigfftmag(:,i)-(del.*alpha.*noiseAvg)),Gain*noiseAvg);
    
    
end

output=signal_recon(sig_recon.^0.5,sigfftphase,win_len,win_overlap);

%diplay computation time
toc;
disp(['Total time for denoising : ' num2str(toc)]);


%de-emphasis
data=filter(1,[1 -0.96],data);
output=filter(1,[1 -0.96],output);


%playing input file
disp('Playing input audio file....');

g=audioplayer(data,Fs);
play(g)
pause(floor(length(data)/Fs));

%playing denoised file 
disp('Playing denoised audio File');
g=audioplayer(output,Fs);
play(g);
pause(floor(length(data)/Fs));


clear all






