%Implement a simple spectral subtraction speech enhancer
%Suppression of Acoustic Noise in Speech Using Spectral Subtraction-S.F Boll
%VAD used is a simple spectral threshold as described in paper.
% need to check performance against LTSV(expected to be better but by how much?)

close all
clear all

tic;
%read noisy data
[data,Fs]=audioread('noisy.wav');

%premphasis
data=filter([1 -0.95],1,data);

%noise silence period
len_ana=0.020; %window length of analysis is 25ms
init_sil=0.25;% initial silence in seconds
overlap_per=0.5; %percentage overlap in windows for enframing 

%Use hanning window as described in paper
win=hamming(len_ana*Fs);
win_len=len_ana*Fs;

%create a function that generates overlapping windows of the input signal
%taking into account the initial silence period
%what do I need?signal, overlap, windowlength,make this a function!!
%need the vector column-wise as fft is computed column-wise by matlab

win_overlap=floor(overlap_per*length(win));

%neglect the last few files less than win
noFr=floor((length(data))/win_overlap-1);

%generate indices of the overlapping matrices
matInd=(repmat(1:length(win),noFr,1)+repmat((0:win_overlap:(noFr-1)*win_overlap)',1,length(win)))';
enframe_sigMat=data(matInd);

clear matInd;

%multiply the frame by the window
enframe_sigMat=enframe_sigMat.*repmat(win,1,noFr);
%clear win;


%now comes the actual processing Shiv !!!- I need to do FFT,compute mag,
% bias, halfwave, reduce residual and them VAD using spectral subtraction
% need to keep updating noise spectral average estimate as I detect using
% VAD

%compute FFT using  the same FFT length as window length
sigfftMat=fft(enframe_sigMat,win_len);
%fft is symmetric so compute magnitude of half the length of each column+DC

sigfftmag = abs(sigfftMat(1:floor(size(sigfftMat,1)/2)+1,:));
sigfftphase= angle(sigfftMat(1:floor(size(sigfftMat,1)/2)+1,:));


%compute the noise average from the initial silence segments
no_iniSilFr=floor(init_sil*Fs/win_overlap-1);

%compute mean rowwise for the intial silence period. This can be updated
%when the VAD detects the period. To keep track of the length of the
%silence period to update it, I use sil_len

noiseAvg=mean(sigfftmag(:,1:no_iniSilFr),2);
noise_len=no_iniSilFr;




%Perform Magnitude Averaging for all your speech frames 
sigfftmag_avg=sigfftmag;
for i=3:noFr-2
   %sigfftmag_avg(:,i)=0.09*sigfftmag(:,i-2)+0.25*sigfftmag(:,i-1)+0.32*sigfftmag(:,i)+0.25*sigfftmag(:,i+1)+0.09*sigfftmag(:,i+2);
   %sigfftmag_avg(:,i)=0.2*sigfftmag(:,i-2)+0.2*sigfftmag(:,i-1)+0.2*sigfftmag(:,i)+0.2*sigfftmag(:,i+1)+0.2*sigfftmag(:,i+2);
     sigfftmag_avg(:,i)=0.25*sigfftmag(:,i-1)+0.5*sigfftmag(:,i)+0.25 *sigfftmag(:,i+1);
end

%Intialize variables for denoising
Gain=0.01;
dbthresh=6;
temp_thresh=8;
counter=9;
isspeech=0;
sig_recon=zeros(floor(win_len/2)+1,noFr);
maxResidue=zeros(floor(win_len/2)+1,1);

%process complete signal frame wise

for i=1:noFr
    %VAD
    [isspeech,counter]=spectralDistVad(sigfftmag(:,i),noiseAvg,counter,dbthresh,temp_thresh);
    if(isspeech==0)
    %update estimate of average noise h
        noise_len=noise_len+1;
        noiseAvg=((noise_len-1)*noiseAvg+sigfftmag(:,i))./noise_len;
    %update maximum noise residue
        maxResidue=max(maxResidue,sigfftmag_avg(:,i)-noiseAvg);
        sig_recon(:,i)=Gain*sigfftmag(:,i);
    
    else
        sig_recon(:,i)=sigfftmag_avg(:,i)-noiseAvg;
        if(i>1 && i<noFr)
            fl=(sig_recon(:,i)< maxResidue);
            sig_recon(fl==1,i)=min([(sigfftmag_avg(fl==1,i-1)-noiseAvg(fl==1)),(sigfftmag_avg(fl==1,i)-noiseAvg(fl==1)),(sigfftmag_avg(fl==1,i+1)-noiseAvg(fl==1))],[],2);
        end
        fl=(sig_recon(:,i)<=0);
        sig_recon(fl,i)=0;
        
    end
    
end


output=OverlapAdd2(sig_recon,sigfftphase,win_len,win_overlap);

%diplay computation time
toc;
disp(['Total time for denoising : ' num2str(toc)]);


%de-emphasis
data=filter(1,[1 -0.95],data);
output=filter(1,[1 -0.95],output);


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






