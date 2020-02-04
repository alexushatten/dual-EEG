function mean_x = alt_mean(x)
% inspired from fastcorr.m by Liber Eleutherios
%
% ... When called only once, it is slower than Matlab's builtin mean.m
% ... When called 1'000'000 times or more, it becomes about 18-20 times 
% faster than mean.m !
%
% Renaud Marquis @FBMlab, June 2018

mean_x = x(1);
for i = 2:length(x)
    delta_x = x(i) - mean_x;
    mean_x = mean_x + delta_x / i; % originally / (i + 1.0) in C but i was equal to 1 (and indexing starts at 0 in C)
end

end