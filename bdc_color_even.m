function c_out = bdc_color_even(z, lim_in, c, c_nr)
    if(nargin < 2)
        lim_in = [ 1/3, 2/3 ];
    end
    if(nargin < 3)
        c = [ 0 0 1; 0 1 0; 1 1 0; 1 0 0; 1 1 1 ];
    end
    if(nargin < 4)
        c_nr = 4000;
    end
    
    nr = numel(lim_in) + 1;
    assert(size(c, 1) >= numel(lim_in) + 2, 'Too many limits');
    assert(numel(c) == 3 * size(c, 1), 'c must be [:] * 3');
    assert(max(max(c)) <= 1 && min(min(c)) >= 0, 'Elements of c must be in [0, 1]');
    
    tmp = sort(reshape(z, numel(z), 1));
    tmp_nr = numel(tmp);
    tmp = (tmp - tmp(1)) * (1 / (tmp(end) - tmp(1)));
    lim = tmp(round(tmp_nr * lim_in));
    lim = reshape(lim, 1, numel(lim_in));

    c_out = zeros(c_nr, 3);
    
    lim = round(c_nr * [0, lim, 1]);
    
    for n = 1:nr
        assert(lim(n) < lim(n + 1), 'lim values must be in increasing order in the range ]0,1[');
        c_out((lim(n) + 1):(lim(n + 1)), :) = bdc_color(lim(n + 1) - lim(n), c(n, :), c(n + 1, :));
    end
end

