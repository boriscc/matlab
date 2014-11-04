function s = get_sample_stat_norm(x)
    s = struct();
    s.n = length(x);
    s.pop_mean = mean(x);
    s.sample_mean = s.pop_mean;
    s.m = s.pop_mean;
    s.m2 = sum((x-s.m).^2) / s.n;
    s.m3 = sum((x-s.m).^3) / s.n;
    s.m4 = sum((x-s.m).^4) / s.n;
    s.k1 = s.m;
    s.k2 = s.n * s.m2 / (s.n-1);
    s.k3 = s.m3 * s.n^2/((s.n-1)*(s.n-2));
    s.k4 = s.n^2 * ((s.n+1)*s.m4 - 3*(s.n-1)*s.m2^2) / ((s.n-1)*(s.n-2)*(s.n-3));
    K = sqrt((s.n-1)/2)*exp(gammaln(0.5*s.n-0.5) - gammaln(0.5*s.n));
    
    s.sample_var = s.m2;
    s.sample_stddev = sqrt(s.sample_var);
    
    s.pop_var = s.k2;
    s.pop_stddev = K * sqrt(s.pop_var);
    
    s.sample_mean_stddev = s.sample_stddev / sqrt(s.n);
    s.sample_mean_var = s.sample_var / s.n;
    
    s.pop_mean_stddev = s.pop_stddev / sqrt(s.n);
    s.pop_mean_var = s.pop_var / s.n;
    
    s.pop_stddev_stddev = s.pop_stddev * K * sqrt(K*K-1);
    s.pop_stddev_var = s.pop_stddev_stddev^2;
    
    s.pop_var_stddev = s.pop_var * sqrt(2 / (s.n-1));
    s.pop_var_var = s.pop_var_stddev^2;
    
    % Sample skewness
    s.sample_skew = s.k3 / s.k2^1.5;
    s.sample_skew_var = 6*s.n*(s.n-1) / ((s.n-2)*(s.n+1)*(s.n+3));
    s.sample_skew_stddev = sqrt(s.sample_skew_var);
    % Population skewness
    s.pop_skew = s.m3 / s.k2^1.5;
    s.pop_skew_var = s.sample_skew_var; % Incorrect, pop_skew_var < sample_skew_var so I use this as approx
    s.pop_skew_stddev = sqrt(s.pop_skew_var);
    % Sample excess kurtosis
    s.sample_exkurt = s.m4 / s.m2^2 - 3;
    s.sample_exkurt_var = 24*s.n*(s.n-1)^2 / ((s.n-3)*(s.n-2)*(s.n+3)*(s.n+5));
    s.sample_exkurt_stddev = sqrt(s.sample_exkurt_var);
    % Unbiased (in the norm case) estimator of the population excess
    % kurtosis
    s.pop_exkurt = s.k4 / s.pop_var^2;
    s.pop_exkurt_var = s.sample_exkurt_var; % Incorrect, I use this as an approximation
    s.pop_exkurt_stddev = sqrt(s.sample_exkurt_var);
end

