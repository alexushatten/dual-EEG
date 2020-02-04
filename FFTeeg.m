function [p,f]=FFTeeg(EEG,fs,epoch,novlap,display)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimates Power Spectral Density by periodogram method using fft()function
%
% Usage:  [p,f]=FFTeeg(EEG,fs,epoch,novlap,display)
%
% INPUTS: 
% -eeg: Time history vector to analyze.
% 
% -fs: Sample frequency of data, used to normalize output
%       and generate frequency vector.
% -epoch: duration of epoch in seconds, will determine nfft that
%       needs to be an even number, but does not
%       have to be a power of two (MATLAB will use
%       a slower DFT routine if not a power of two,
%       see 'help fft' for details).
%
% -novlap: percentage overlap, for example nfft=1024 with novlap=512 is
%       a 50% overlap.
% - notch: if the signal was notch filtered
%  -fres: frequ resolution, smoothes curve by grouping frequ bins
%
% OUTPUTS 
% -p: Power spectral density in units of [y units]
%       squared per [fs units], for example g^2/Hz.
% -f: Frequency bins.
%
% PROCESSING
%   - wndowing: Hanning window of length nfft to each epoch 
%   - oarms: Overall rms value, square root of area under f-p
%       curve.
%
% DEVELOPMENT
%   - Michel Speiser 2010, lasr package
%   - Max Baud July 2015, same concept, rewrote to have epoch by epoch
%   power
%
% TO DOs
%   - [] implement frequ resolution option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs
if numel(EEG)>length(EEG)
    error('Input only one channel at a time')
end
if size(EEG,1)<size(EEG,2)
    EEG=EEG'; % column vector
end


% time series
nfft=floor(epoch*fs/2)*2; % epoch*fs, make sure it is an even integer
novlap=floor(novlap*nfft); % overlap in datapoints
nepoch = floor((length(EEG)-nfft)/(nfft-novlap))+1 ;% number of epochs to loop through

% frequency bins
df = fs/nfft ;
f = (0:df:fs/2)';
        
% Hanning window
w=hann(nfft); 
w_sc = sum(w.^2) ; %normalization factor for window
% use  w_sc = nfft ; if no windowing

% Initialize 
        n1 = 1 ;
        n2 = nfft ;
        dn = nfft-novlap ;
        p=zeros(nepoch,nfft);
        arg_sum = zeros([nfft 1]) ; % Initialize psd summation storage variable
        
        
%%%% Main loop
for e=1:nepoch
    e_eeg = EEG(n1:n2) ; % Extract current epoch eeg
    e_eeg = e_eeg - mean(e_eeg) ; % Centered 
    e_eeg = e_eeg.*w; % Apply window
    
    e_fft = fft(e_eeg,nfft) ; % FFT of epoch
    e_fft = abs(e_fft).^2 ; % Modulus squared of FFT
    p(e,:)=e_fft;
    % Increment epoch index
    n1 = n1+dn ; 
    n2 = n2+dn ;
end
%%%% End of main loop

% Normalize
% Power spectrum is symmetric about Nyquist frequency, use lower half
%and multiply by 4. Then divide by 2 to convert peak^2 to rms^2. Thus
%scale factor is 2.
p = 2*p(:,1:nfft/2+1) ;
p = p/w_sc ;% correct for window power
p = p/fs ; % Normalize spectral density for sampling units
%oarms = sqrt(sum(p.*df)) ; % Calculate overall rms level from area under PSD curve


%%% FIGURE   
        if display
            if e>1
                p2=nanmean(p);
                col=cbrewer('seq','Blues',size(p,1));
            else col=[1 0 0; 0 0 1];
            end
            figure
            set(gca,'colororder',col,'yscale','log')
            hold on
            plot(f,p)
            if e>1;plot(f,p2,'r-.','linewidth',2);end
            
            xlim([0 max(f)])
%             ylim([0 1.2*max(p(:))])
            xlabel('Frequency (Hz)')
            ylabel('Mean Power (µV^2/frequ)')
            
        end
end