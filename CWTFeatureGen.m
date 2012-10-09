function [out] = CWTFeatureGen(myCWT)
% [] = CWTFeatureGen()
% Input:
%
%
% Output:
% out.cwt(i).features =
%
% Description:
% Can the information from a CWT be transformed into a set of features that
% can be used to distinguish singals from one another? This function will
% create the features that will be used to answer this question.
%
% Other Notes:
%
out = myCWT;
for i=1:length(myCWT.cwt)
    low = 1:round(length(myCWT.scales)/3);
    med = (low(end)+1):round(length(myCWT.scales)*2/3);
    high = (med(end)+1):length(myCWT.scales);
    peakinfo = [myCWT.cwt(i).ridgepeaks.scaleindex, myCWT.cwt(i).ridgepeaks.time, myCWT.cwt(i).ridgepeaks.cfs];
    if isempty(peakinfo) %no peaks were found...
        low.num = 0;
        med.num = 0;
        high.num = 0;
        low.meancfs = 0;
        med.meancfs = 0;
        high.meancfs = 0;
        out.cwt(i).features.low = low;
        out.cwt(i).features.med = med;
        out.cwt(i).features.high = high;
        continue;
    end
    peakinfo = sortrows(peakinfo,1);
    lowindtemp = find(peakinfo(:,1)<=low(end),1,'last');
    medindtemp = find((peakinfo(:,1)>low(end))&(peakinfo(:,1)<=med(end)),1,'last');
    lowind = zeros(size(myCWT.cwt(i).ridgepeaks.scaleindex));
    if ~isempty(lowindtemp)
        lowind(1:lowindtemp) = 1;
    else
        lowindtemp = 0;
    end
    medind = zeros(size(myCWT.cwt(i).ridgepeaks.scaleindex));
    if ~isempty(medindtemp)
        medind((lowindtemp+1):medindtemp) = 1;
    else
        medindtemp = lowindtemp;
    end
    highind = zeros(size(myCWT.cwt(i).ridgepeaks.scaleindex));
    if (medindtemp+1)<=length(highind)
        highind((medindtemp+1):end) = 1;
    end
    
    lowind = logical(lowind);
    medind = logical(medind);
    highind = logical(highind);
    %count the number of peaks
    low.num = sum(lowind);
    med.num = sum(medind);
    high.num = sum(highind);
    %calculate the average value of these peaks
    cfsarray = peakinfo(:,3);
    if low.num ~= 0
        low.meancfs = mean(cfsarray(lowind));
    else
        low.meancfs = 0;
    end
    if med.num ~=0
        med.meancfs = mean(cfsarray(medind));
    else
        med.meancfs = 0;
    end
    if high.num ~= 0
        high.meancfs = mean(cfsarray(highind));
    else
        high.meancfs = 0;
    end
    out.cwt(i).features.low = low;
    out.cwt(i).features.med = med;
    out.cwt(i).features.high = high;
end

out.vector = [out.cwt(1).features.low.num, out.cwt(1).features.low.meancfs, out.cwt(1).features.med.num, out.cwt(1).features.med.meancfs, out.cwt(1).features.high.num, out.cwt(1).features.high.meancfs, ...
    out.cwt(2).features.low.num, out.cwt(2).features.low.meancfs, out.cwt(2).features.med.num, out.cwt(2).features.med.meancfs, out.cwt(2).features.high.num, out.cwt(2).features.high.meancfs, ...
    out.cwt(3).features.low.num, out.cwt(3).features.low.meancfs, out.cwt(3).features.med.num, out.cwt(3).features.med.meancfs, out.cwt(3).features.high.num, out.cwt(3).features.high.meancfs];