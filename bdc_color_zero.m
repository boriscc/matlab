function col = bdc_color_zero(Z, col_nr, black_nr)
    Z_min = min(min(Z));
    Z_max = max(max(Z));
    if nargin < 2
        col_nr = 2048;
    end
    if nargin < 3
        black_nr = 16;
    end
    col = zeros(col_nr, 3);
    negneg = [0 0 1];
    neg = [0 1 0];
    zero = [0 0 0] + 0.0;
    pos = [1 1 0];
    pospos = [1 0 0];
    mode = 'linear';
    power = 1;
    if Z_min < 0 && Z_max > 0
        nr_neg = round(-col_nr * Z_min / (Z_max - Z_min) - black_nr);
        nr_pos = col_nr - nr_neg - black_nr * 2;
        if nr_neg > 0
            col(1:nr_neg, :) = bdc_color(nr_neg, negneg, neg, power, 1, mode);
        end
        if nr_neg + black_nr > 0
            if nr_neg < 0
                start = 0;
            else
                start = nr_neg;
            end
            col((start + 1):(start + black_nr), :) = bdc_color(black_nr, neg, zero, power, 1, mode);
        end
        if nr_neg + black_nr < col_nr
            idx = (nr_neg + black_nr + 1):min(nr_neg + 2 * black_nr, col_nr);
            col(idx, :) = bdc_color(length(idx), zero, pos, power, 1, mode);
            if nr_pos > 0
                col((end - nr_pos + 1):end, :) = bdc_color(nr_pos, pos, pospos, power, 1, mode);
            end
        end
    elseif Z_min < 0 && Z_max <= 0
        if -Z_min * 0.01 > -Z_max
            col(1:(col_nr - black_nr), :) = bdc_color(col_nr - black_nr, negneg, neg, power, 1, mode);
            col((col_nr - black_nr + 1):end, :) = bdc_color(black_nr, neg, zero, power, 1, mode);
        else
            col(1:end, :) = bd_color(col_nr, negneg, neg, power, 1, mode);
        end
    elseif Z_min >= 0 && Z_max > 0
        if Z_max * 0.01 > Z_min
            col(1:black_nr, :) = bdc_color(black_nr, zero, pos, power, 1, mode);
            col((black_nr + 1):end, :) = bdc_color(col_nr - black_nr, pos, pospos, power, 1, mode);
        else
            col(1:end, :) = bdc_color(col_nr, pos, pospos, power, 1, mode);
        end
    else
        error('Not implemented');
    end
end
