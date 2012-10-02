function [out] = ridgefinder(myCWT)
% [] = CWTmodded()
% Input:
% myCWT = the cwt struct output by CWTmodded().
% myCWT.scales = the wavelet scales used in the transform
% myCWT.cwt
% myCWT.cwt.cfs = the wavelet coefficients from the transform
% myCWT.cwt.type = a string with the type of wavelet used in the transform
% myCWT.cwt(i) = there will be 3 different wavelet transforms computed: mexican-hat, derivative-of-gaussian, and the real-valued Morlet wavelet.
%
% Output:
% out.rigdgmap = a mapping of a ridge label onto the scalogram
% out.ridge = 
%
% Description:
% The signal that has been transformed into a wavelet domain can be thought
% of as a surface with peaks and valleys. The peaks of a signal often
% persists across several scales, creating ridges. These ridges can be
% identified and the peak of each ridge, too. These ridges and peaks offer
% a more condensed description of a wavelet transform and can indicate the
% presence of an interesting feature in the original waveform.
%
%
% Other Notes:
%
% 1. Find the peaks within each scale
out = myCWT;
numScales = length(myCWT.scales);
numWavelets = length(myCWT.cwt);
numTime = size(myCWT.cwt(1).cfs,2);
wavelet_peaks = cell(1,numWavelets);
for j = 1:numWavelets
    wavelet_peaks{j} = cell(1,numScales);
    for k = 1:numScales
        wavelet_peaks{j}{k} = first_pass_peak_detection(myCWT.cwt(j).cfs(k,:), out.scales(k)*2+1, 0);
    end
end
%3b.ii. Find the ridges
wavelet_peaks_alias = cell(1,numWavelets);
gap_limit = 2; %connect ridges that are even separated by "gap_limit-1" scales;
%%%%% ridge finding code
for z = 1:numWavelets
    ridge_map = zeros(numScales,numTime);
    wavelet_peaks_alias{z} = cell(1,numScales);
    for i=1:numScales
        wavelet_peaks_alias{z}{i} = zeros(size(wavelet_peaks{z}{i}));
    end
    for i=1:length(wavelet_peaks{z}{end})
        ridge_map(end,wavelet_peaks{z}{end}(i)) = i;
        wavelet_peaks_alias{z}{end}(i) = i;
    end
    ridge_counter = length(wavelet_peaks{z}{end}); %keeps track of the total number of ridges
    for i=numScales:-1:(gap_limit+1)
        for j=1:length(wavelet_peaks{z}{i})
            for h=1:gap_limit
                %Search for peaks within the window size for scale i.
                low_bnd = wavelet_peaks{z}{i}(j) - 2*myCWT.scales(i-h);
                up_bnd = wavelet_peaks{z}{i}(j) + 2*myCWT.scales(i-h);
                low_set = wavelet_peaks{z}{i-h}>low_bnd;
                up_set = wavelet_peaks{z}{i-h}<up_bnd;
                %If a peak is found add it to the growing ridge
                if any(low_set.*up_set)
                    ridge_set = low_set.*up_set;
                    if sum(ridge_set)==1
                        for k=1:length(ridge_set)
                            if ridge_set(k) && (wavelet_peaks_alias{z}{i-h}(k)==0 && ...
                                    sum(wavelet_peaks_alias{z}{i-h}==wavelet_peaks_alias{z}{i}(j))==0)
                                ridge_map(i-h,wavelet_peaks{z}{i-h}(k)) = wavelet_peaks_alias{z}{i}(j);
                                wavelet_peaks_alias{z}{i-h}(k) = wavelet_peaks_alias{z}{i}(j);
                            end
                        end
                    else
                        my_min = inf;
                        for k=1:length(ridge_set)
                            if ridge_set(k) && wavelet_peaks_alias{z}{i-h}(k)==0
                                temp_ind=wavelet_peaks_alias{z}{i-(h-1)}==wavelet_peaks_alias{z}{i}(j);
                                penultimate_ridge_position=wavelet_peaks{z}{i-(h-1)}(temp_ind);
                                temp_min = abs(wavelet_peaks{z}{i-h}(k)-penultimate_ridge_position);
                                if temp_min<my_min
                                    k_min = k;
                                    my_min = temp_min;
                                end
                            end
                        end
                        if (sum(wavelet_peaks_alias{z}{i-h}==wavelet_peaks_alias{z}{i}(j))==0) && (k_min<=k)
                            ridge_map(i-h,wavelet_peaks{z}{i-h}(k_min)) = wavelet_peaks_alias{z}{i}(j);
                            wavelet_peaks_alias{z}{i-h}(k_min) = wavelet_peaks_alias{z}{i}(j);
                        end
                    end
                end
            end
        end
        %Start new ridges for the peaks that were not assigned to ridges in the
        %scale below.
        new_ridge_set = wavelet_peaks_alias{z}{i-1}==0;
        for j=1:length(wavelet_peaks_alias{z}{i-1})
            if new_ridge_set(j)
                ridge_counter = ridge_counter + 1;
                wavelet_peaks_alias{z}{i-1}(j) = ridge_counter;
                ridge_map(i-1,wavelet_peaks{z}{i-1}(j)) = ridge_counter;
            end
        end
    end
    out.cwt(z).ridgemap = ridge_map;
    out.cwt(z).ridge = ridgemap2ridge(ridge_map);
end

function [out]=first_pass_peak_detection(x,winw_size,s)
%There are many ways to find peaks. The wavelet method is powerful at
%detecting peaks but is actually dependent on simpler peak detection
%methods.

%Simple Peak Finding Method: Scan a signal with a window. Find the max. If the max is greater
%than the left and right endpoints of the window it is a peak candidate.
%The window is centered at each point of a waveform, so a wider peak
%candidate will recieve more votes in a sense. If a window has more than
%one point with the max value then the left most index is used. A downside
%to this method is that it is sensitive to the size of the window. However,
%this seeming limitation is actually put to use in the wavelet method by
%scaling the window along with the wavelet.
length_x = length(x);
if size(x,2) == 1
    x=x';
end
if winw_size<3
    winw_size = 3;
    warning('win_size:small','The window size for first pass peak detection may be too small and therefore yield illogical results.');
end
pad=(winw_size-1)/2;
[x_padded,window_padded] = deal(zeros(1,length_x+winw_size-1));
window_padded(1:winw_size) = ones(1,winw_size);
x_padded(pad+1:end-pad) = x;
out=zeros(size(x));
for i=1:length_x
    seg = x_padded(logical(window_padded));
    [~,ind] = max(seg);
    if seg(ind) > seg(1) && seg(ind) > seg(end)
        out(i) = i + ind - pad - 1;
    end
    window_padded = circshift(window_padded,[0, 1]);
end
out(out==0) = [];
peak_candidates = unique(out);
peak_elected = peak_candidates; %This vector will be trimmed below
%Tally votes; if a peak candidate has less than or equal to pad votes they
%are disqualified
for i=length(peak_candidates):-1:1
    temp = out==peak_candidates(i);
    if sum(temp)<=pad/2
        peak_elected(i) = [];
    end
end
%Weed out unqualified peaks and peaks of questionable nature.
%If a peak has a negative wavelet coefficient value it is disqualified and
%if a peak is the first or last point of data it is disqualified.
for i=length(peak_elected):-1:1
    if x(peak_elected(i))<0
        peak_elected(i) = [];
    elseif peak_elected(i) == 1
        peak_elected(i) = [];
    elseif peak_elected(i) == length_x
        peak_elected(i) = [];
    end
end

for i=length(peak_elected):-1:1
    if x(peak_elected(i))<s
        peak_elected(i) = [];
    end
end
out = peak_elected;

function [ridge] = ridgemap2ridge(rmap)
numRidge = max(max(rmap));
ridge = cell(numRidge,1);
for i = 1:numRidge
    ridge{i} = find(rmap==i);
end

function [peaks] = ridgepeaks(cwt)
numRidge = length(cwt.ridge);
for i=1:numRidge
    ridgecfs = cwt.cfs(cwt.ridge{i});
    peak_index = watershed(ridgecfs); %easy way to find peaks
    peak_index = find(~peak_index);
end