function [out] = CWTmodded(mySignal)
% [] = CWTmodded()
% Input:
% mySignal = the 1-D signal for which to calculate a transform.
%
% Output:
% out = a struct that contains the wavelet transform(s) and other information.
% out.scales = the wavelet scales used in the transform
% out.cwt
% out.cwt.cfs = the wavelet coefficients from the transform
% out.cwt.type = a string with the type of wavelet used in the transform
% out.cwt(i) = there will be 3 different wavelet transforms computed: mexican-hat, derivative-of-gaussian, and the real-valued Morlet wavelet.
%
% Description:
% A wavelet transform can be calculated using an FFT algorithm in MATLAB
% using the cwtft(). However, the scales that are returned must be evenly
% spaced. It is sometimes convenient to return scales that are spaced more
% exponentially, so this function will do that automatically based upon the
% length of the signal.
%
% In addition, this function applies a tukey window to the input signal in
% an effort to suppress edge effects -- artificial sharp transitions at
% the beginning and end of a signal -- of the wavelet transform.
%
% Other Notes:
% Choosing the right scales to investigate can be a challenge, because it
% can feel subjective. One way is to include every integer scale up to the
% length of the signal, but this is probably too much information. Another
% is to choose scales on an exponential/log scale, but this might skip over
% some interesting scales. We'll use a hybrid between the two: filling
% in an exponential scale with uniform spacing. Hopefully this compromise
% will deliver detail across several orders of magnitude.

% 1. The wavelet transform can be computed using a FFT algorithm. In order to leverage the advantages of the FFT (speed), expand the number of data points to the nearest power of two of 150% of the original signal length. The 150% is needed for effective windowing.
numTime = length(mySignal);
fftLen= 2^(ceil(log(numTime*2)/log(2)));
timeDiff = fftLen - numTime;
if mod(timeDiff,2)
    %is odd
    signalPow2 = padarray(mySignal,[0 (timeDiff-1)/2],'replicate');
    signalPow2(end+1) = signalPow2(end);
    originalSignalIndLeft = 1+(timeDiff-1)/2;
else
    %is even
    signalPow2 = padarray(mySignal,[0 timeDiff/2],'replicate');
    originalSignalIndLeft = 1+(timeDiff)/2;
end
% 2. Apply a Tukey window, a tapered cosine, to suppress edge effects. The expanded signal will be lowered to zero while the original signal remains the same.
percentSignalExpansion = timeDiff/fftLen;
tukeyWindow = tukeywin(fftLen, percentSignalExpansion);
signalWin = signalPow2.*tukeyWindow';

% 3. Determine the scales that will have wavelet coefficients.
halfSupport = numTime/8; %Based upon the effective support of a wavelet. See 'waveinfo()'.
pow102 = ceil(log(halfSupport/10)/log(2)); % i.e. 1,2,3...,8,9,10,12,14,16...,26,28,30,34,38,42...,62,66,70...
if pow102<2
    pow102=2;
end
wavelet_scales = cell(1,pow102);
for i=1:pow102
    wavelet_scales{i} = (1:10)*2^(i-1)+10*(2^(i-1)-1);
end
if(wavelet_scales{pow102}(1)>halfSupport)
    wavelet_scales{pow102}(2:end) = [];
else
    wavelet_scales{pow102}(wavelet_scales{pow102}>halfSupport) = [];
end
out.scales = cell2mat(wavelet_scales);
numScales = length(out.scales);

% 4a. Calculate the mexican hat wavelet transform
out.cwt(1).type = 'mexh';
cwtftOutputTemp = cwtft(signalWin,'scales',wavelet_scales{1},'wavelet','mexh');
cwtftOutputTemp2 = zeros(numScales,fftLen);
cwtftOutputTemp2(1:10,:) = cwtftOutputTemp.cfs;
if (pow102>2)
    for i=2:(pow102-1)
        cwtftOutputTemp = cwtft(signalWin,'scales',wavelet_scales{i},'wavelet','mexh');
        cwtftOutputTemp = cwtftOutputTemp.cfs;
        cwtftOutputTemp2((pow102-1)*10+1:(pow102-1)*10+10,:) = cwtftOutputTemp;
    end
end
cwtftOutputTemp = cwtft(signalWin,'scales',wavelet_scales{pow102},'wavelet','mexh');
cwtftOutputTemp = cwtftOutputTemp.cfs;
cwtftOutputTemp2((end-length(wavelet_scales{pow102})+1):end,:) = cwtftOutputTemp;
cwtftOutputTemp2 = real(cwtftOutputTemp2);
for i = 1:numScales
    cwtftOutputTemp2(i,:) = cwtftOutputTemp2(i,:)/sqrt(out.scales(i));
end
out.cwt(1).cfs = cwtftOutputTemp2(:,originalSignalIndLeft:(originalSignalIndLeft+numTime-1));

% 4b. Calculate the derivate of gaussian wavelet transform
out.cwt(2).type = 'dog';
cwtftOutputTemp = cwtft(signalWin,'scales',wavelet_scales{1},'wavelet',{'dog',1});
cwtftOutputTemp2 = zeros(numScales,fftLen);
cwtftOutputTemp2(1:10,:) = cwtftOutputTemp.cfs;
if (pow102>2)
    for i=2:(pow102-1)
        cwtftOutputTemp = cwtft(signalWin,'scales',wavelet_scales{i},'wavelet',{'dog',1});
        cwtftOutputTemp = cwtftOutputTemp.cfs;
        cwtftOutputTemp2((pow102-1)*10+1:(pow102-1)*10+10,:) = cwtftOutputTemp;
    end
end
cwtftOutputTemp = cwtft(signalWin,'scales',wavelet_scales{pow102},'wavelet',{'dog',1});
cwtftOutputTemp = cwtftOutputTemp.cfs;
cwtftOutputTemp2((end-length(wavelet_scales{pow102})+1):end,:) = cwtftOutputTemp;
cwtftOutputTemp2 = real(cwtftOutputTemp2);
for i = 1:numScales
    cwtftOutputTemp2(i,:) = cwtftOutputTemp2(i,:)/sqrt(out.scales(i));
end
out.cwt(2).cfs = cwtftOutputTemp2(:,originalSignalIndLeft:(originalSignalIndLeft+numTime-1));

% 4c. Calculate the real-valued morlet wavelet transform
out.cwt(3).type = 'morlex';
cwtftOutputTemp = cwtft(signalWin,'scales',wavelet_scales{1},'wavelet','morlex');
cwtftOutputTemp2 = zeros(numScales,fftLen);
cwtftOutputTemp2(1:10,:) = cwtftOutputTemp.cfs;
if (pow102>2)
    for i=2:(pow102-1)
        cwtftOutputTemp = cwtft(signalWin,'scales',wavelet_scales{i},'wavelet','morlex');
        cwtftOutputTemp = cwtftOutputTemp.cfs;
        cwtftOutputTemp2((pow102-1)*10+1:(pow102-1)*10+10,:) = cwtftOutputTemp;
    end
end
cwtftOutputTemp = cwtft(signalWin,'scales',wavelet_scales{pow102},'wavelet','morlex');
cwtftOutputTemp = cwtftOutputTemp.cfs;
cwtftOutputTemp2((end-length(wavelet_scales{pow102})+1):end,:) = cwtftOutputTemp;
cwtftOutputTemp2 = real(cwtftOutputTemp2);
for i = 1:numScales
    cwtftOutputTemp2(i,:) = cwtftOutputTemp2(i,:)/sqrt(out.scales(i));
end
out.cwt(3).cfs = cwtftOutputTemp2(:,originalSignalIndLeft:(originalSignalIndLeft+numTime-1));