function [ bin_peak_max, val_peak_max, peak_val, peak_abs ] = peak_detector_general( signal, method, opt_arg )
%peak_detector_general detects peaks within the signal acording to the
%specified method.
% signal must be a vector
% method must be a string
%   Possible methods are :
%                 - mean : Threshold of mean(signal) + x*ste(signal)
%                 (opt_arg : x)
%                 - thr  : Threshold the signal by thr value. If thr < 0,
%                 peaks < thr are considered. opt_arg = thr.
%                 
%
% opt_arg is specific to the method used
%
% Output :
%   - bin_peak_max : bin of the maximal value within each peak
%   - val_peak_max : maximum of each peak
%   - peak_val : values of each peak
%   - peak_abs : first and last bin of each peak


switch method
    case 'mean'
        if ~exist('opt_arg','var'),
            opt_arg = 4;
        end
        threshold = mean(signal) + opt_arg*std(signal);
        filt_sign = signal >= threshold;
        dfilt_sign = diff(filt_sign);
        
        ind1 = find(dfilt_sign==1);
        
        % Case no peaks
        if isempty(ind1),
            bin_peak_max = [];
            val_peak_max = [];
            peak_val = {};
            peak_abs = [];
            return;
        end
        
        % First peak
        indm1 = find(dfilt_sign==-1,1,'first');
        if ind1(1)>indm1,
            peak_abs(1,1:2) = [1 indm1];
            peak_val{1} = signal(peak_abs(1,1):peak_abs(1,2));
            [val_peak_max(1) bin_peak_max(1)] = max(peak_val{1});
            c = 2;
        else
            c = 1;
        end
        
        % peaks
        for i = 1:length(ind1)-1,
            indm1 = find(dfilt_sign(ind1(i):end)==-1,1,'first') + ind1(i);
            peak_abs(c,1:2) = [ind1(i)+1 indm1-1];
            peak_val{c} = signal(peak_abs(c,1):peak_abs(c,2));
            [val_peak_max(c) bin_peak_max(c)] = max(peak_val{c});
            bin_peak_max(c) = bin_peak_max(c) + ind1(i);
            c=c+1;
        end
        
        % Case one and only peak in the middle
        if isempty(i);
            i=1;
        end
        
        % Last peak
        indm1 = find(dfilt_sign(ind1(end):end)==-1,1,'first') + ind1(i);
        if ~isempty(indm1),
            peak_abs(c,1:2) = [ind1(i)+1 indm1-1];
            peak_val{c} = signal(peak_abs(c,1):peak_abs(c,2));
            [val_peak_max(c) bin_peak_max(c)] = max(peak_val{c});
            bin_peak_max(c) = bin_peak_max(c) + ind1(i);
        else
            indm1 = length(signal);
            peak_abs(c,1:2) = [ind1(i)+1 indm1];
            peak_val{c} = signal(peak_abs(c,1):peak_abs(c,2));
            [val_peak_max(c) bin_peak_max(c)] = max(peak_val{c});
            bin_peak_max(c) = bin_peak_max(c) + ind1(i);
        end
        
        case 'thr'
        
        threshold = opt_arg;
        flagNeg = 0;
        if threshold < 0,
            flagNeg = 1;
            signal = -1 * signal;
            threshold = abs(threshold);
        end
        
        filt_sign = signal >= threshold;
        dfilt_sign = diff(filt_sign);
        
        ind1 = find(dfilt_sign==1);
        ind2 = find(dfilt_sign==-1);
        
        % Case no peaks
        if isempty(ind1),
            bin_peak_max = [];
            val_peak_max = [];
            peak_val = {};
            peak_abs = [];
            return;
        end
        
        flag1 = 0; % Flag signal starting above threshold
        flag2 = 0; % Flag signal ending above threshold
        if length(ind1) == length(ind2) && ind2(1)-ind1(1)<0
            flag1 = 1;
            flag2 = 1;
            indm1 = ind2(1);
            ind2 = ind2(2:end);
            indm2 = ind1(end);
            ind1 = ind1(1:end-1);
        elseif length(ind1)>length(ind2)
            flag2 = 1;
            indm2 = ind1(end);
            ind1 = ind1(1:end-1);
        elseif length(ind1)<length(ind2)
            flag1 = 1;
            indm1 = ind2(1);
            ind2 = ind2(2:end);
        end
        
        % Preallocate
        peak_abs = zeros(length(ind1)+flag1+flag2,2);
        peak_val = cell(length(ind1)+flag1+flag2,1);
        val_peak_max = zeros(length(ind1)+flag1+flag2,1);
        bin_peak_max = zeros(length(ind1)+flag1+flag2,1);
        
        c = 1;
        % First peak
        if flag1,
            peak_abs(1,1:2) = [1 indm1];
            peak_val{1} = signal(peak_abs(1,1):peak_abs(1,2));
            [val_peak_max(1), bin_peak_max(1)] = max(peak_val{1});
            c = c+1;
        end
        
        % peaks
        for i = 1:length(ind1),
            peak_abs(c,1:2) = [ind1(i)+1 ind2(i)];
            peak_val{c} = signal(peak_abs(c,1):peak_abs(c,2));
            [val_peak_max(c), bin_peak_max(c)] = max(peak_val{c});
            bin_peak_max(c) = bin_peak_max(c) + ind1(i);
            c=c+1;
        end 
        
        % Last peak
        if flag2,
            peak_abs(c,1:2) = [indm2+1 length(signal)];
            peak_val{c} = signal(peak_abs(c,1):peak_abs(c,2));
            [val_peak_max(c), bin_peak_max(c)] = max(peak_val{c});
            bin_peak_max(c) = bin_peak_max(c) + indm2;
        end
        
        if flagNeg, % If trig < 0, put val back to normal
            val_peak_max = -1 * val_peak_max;
            peak_val = cellfun(@(x)(-1 .* x),peak_val,'UniformOutput',false);
        end
        
    otherwise
        error('Unknown method');
end

end

