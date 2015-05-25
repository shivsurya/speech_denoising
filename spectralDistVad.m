function [isspeech,counter]=spectralDistVad(sig,noise,counter,dbthresh,temporal_thresh)

%calculate spectral distance
spectral_dist= mean(abs(20*(log10(sig)-log10(noise))));

if (spectral_dist < dbthresh) 
     
    counter=counter+1;
else
    counter=0;
end

if (counter > temporal_thresh) 
    isspeech=0;    
else 
    isspeech=1; 
end

end