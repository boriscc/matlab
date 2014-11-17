function str = get_num_with_err(num, err, significant)
    if nargin < 3
        significant = 2;
    end
    if imag(err) ~= 0
        error('Got imaginary error');
    end
    
    expo = log10(10^significant / err);
    if expo == floor(expo)
        expo = floor(expo) - 1;
    else
        expo = floor(expo);
    end
    err_disp = round(err*10^expo);
    num_disp = round(num*10^expo)*10^-expo;
    if expo >= 0
        form = sprintf('%%+.%df(%%d)', expo);
        str = sprintf(form, num_disp, err_disp);
    else
        str = sprintf('%+*.0f(%d)', significant - expo, num_disp, err_disp*10^-expo);
    end
end
