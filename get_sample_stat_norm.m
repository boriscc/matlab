function s = get_sample_stat_norm(x)
    s = struct();
    s.n = length(x);
    s.m = mean(x);
    s.m2 = sum((x-s.m).^2) / s.n;
    s.m3 = sum((x-s.m).^3) / s.n;
    s.m4 = sum((x-s.m).^4) / s.n;
    s.k1 = s.m;
    s.k2 = s.n * s.m2 / (s.n-1);
    s.k3 = s.m3 * s.n^2/((s.n-1)*(s.n-2));
    s.k4 = s.n^2 * ((s.n+1)*s.m4 - 3*(s.n-1)*s.m2^2) / ((s.n-1)*(s.n-2)*(s.n-3));
    K = sqrt((s.n-1)/2)*exp(gammaln(0.5*s.n-0.5) - gammaln(0.5*s.n));
    s.var = s.k2;
    s.stddev = K * sqrt(s.var);
    s.m_stddev = s.stddev / sqrt(s.n);
    s.m_var = s.var / s.n;
    s.stddev_stddev = s.stddev * K * sqrt(K*K-1);
    s.stddev_var = s.stddev_stddev^2;
    s.var_stddev = s.var * sqrt(2 / (s.n-1));
    s.var_var = s.var_stddev^2;
    s.g1 = s.k3 / s.var^1.5; % Should it be s.k3 / (s.var * s.stddev) ??
    s.g1_var = 6*s.n*(s.n-1) / ((s.n-2)*(s.n+1)*(s.n+3));
    s.g1_stddev = sqrt(s.g1_var);
    s.g2 = s.k4 / s.var^2;
    s.g2_var = 24*s.n*(s.n-1)^2 / ((s.n-3)*(s.n-2)*(s.n+3)*(s.n+5));
    s.g2_stddev = sqrt(s.g2_var);
end

