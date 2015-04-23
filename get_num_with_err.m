function str = get_num_with_err(num, err_down, err_up, num_significant, err_significant)
    sym = 0;
    if nargin < 3
        err_up = err_down;
        err_down = -err_down;
        sym = 1;
        num_significant = Inf;
        err_significant = 2;
    end
    if nargin < 4
        num_significant = Inf;
        err_significant = 2;
    end
    if nargin < 5
        err_significant = 2;
    end
    if imag(err_down) ~= 0 || imag(err_up) ~= 0
        error('Got imaginary error');
    end
    if -err_down == err_up
        sym = 1;
    end
    
    if sym == 0
        assert(err_down <= 0 && err_up >= 0, 'Bad error limits');
    end
    
    significant = min([get_level(max([-err_down, err_up])) - get_level(num) + num_significant, err_significant]);
    
    if significant <= -1
        level = get_level(num);
        num = 10^(level+1-num_significant)*round(num*10^(num_significant-level-1));
        if level+1-num_significant >= 0
            str = sprintf('%+.0f', num);
        else
            form = sprintf('%%+.%df', num_significant - level - 1);
            str = sprintf(form, num);
        end
        return;
    end
    
    if sym == 1
        err = err_up;
        assert(err > 0, 'Bad error');

        [ expo, err_disp, num_disp ] = get_disp(significant, err, num);
        if get_level(err_disp) == err_significant
            significant = significant - 1;
            [ expo, err_disp, num_disp ] = get_disp(significant, err, num);
        end
        if expo >= 0
            if err_disp == 0
                form = sprintf('%%+.%df', expo);
                str = sprintf(form, num_disp);
            else
                form = sprintf('%%+.%df(%%d)', expo);
                str = sprintf(form, num_disp, err_disp);
            end
        else
            if err_disp == 0
                str = sprintf('%+*.0f', significant - expo, num_disp);
            else
                str = sprintf('%+*.0f(%d)', significant - expo, num_disp, err_disp*10^-expo);
            end
        end
    else
        assert(err_up > err_down, 'Must have one error with |err| > 0');
        if err_up == 0
            [ expo, err_disp, num_disp ] = get_disp(significant, -err_down, num);
            if get_level(err_disp) == err_significant
                significant = significant - 1;
                [ expo, err_disp, num_disp ] = get_disp(significant, -err_down, num);
            end
            if expo >= 0
                if err_disp == 0
                    form = sprintf('%%+.%df', expo);
                    str = sprintf(form, num_disp);
                else
                    form = sprintf('%%+.%df_{(-%%d)}^{(+0)}', expo);
                    str = sprintf(form, num_disp, err_disp);
                end
            else
                if err_disp == 0
                    str = sprintf('%+*.0f', significant - expo, num_disp);
                else
                    str = sprintf('%+*.0f_{(-%d)}^{(+0)}', significant - expo, num_disp, err_disp*10^-expo);
                end
            end
        elseif err_down == 0
            [ expo, err_disp, num_disp ] = get_disp(significant, err_up, num);
            if get_level(err_disp) == err_significant
                significant = significant - 1;
                [ expo, err_disp, num_disp ] = get_disp(significant, err_up, num);
            end
            if expo >= 0
                if err_disp == 0
                    form = sprintf('%%+.%df', expo);
                    str = sprintf(form, num_disp);
                else
                    form = sprintf('%%+.%df_{(-0)}^{(+%%d)}', expo);
                    str = sprintf(form, num_disp, err_disp);
                end
            else
                if err_disp == 0
                    str = sprintf('%+*.0f', significant - expo, num_disp);
                else
                    str = sprintf('%+*.0f_{(-0)}^{(+%d)}', significant - expo, num_disp, err_disp*10^-expo);
                end
            end
        else
            if -err_down < err_up
                err_max = err_up;
            else
                err_max = -err_down;
            end
            [ expo, err_disp, num_disp ] = get_disp(significant, err_max, num);
            if get_level(err_disp) == err_significant
                significant = significant - 1;
                [ expo, err_disp, num_disp ] = get_disp(significant, err_max, num);
            end
            if -err_down > err_up
                err_down_disp = err_disp;
                err_up_disp = round(err_up*10^expo);
            else
                err_up_disp = err_disp;
                err_down_disp = -round(err_down*10^expo);
            end
            if expo >= 0
                if err_up_disp == 0 && err_down_disp == 0
                    form = sprintf('%%+.%df', expo);
                    str = sprintf(form, num_disp);
                elseif err_up_disp == err_down_disp
                    form = sprintf('%%+.%df(%%d)', expo);
                    str = sprintf(form, num_disp, err_up_disp);
                else
                    form = sprintf('%%+.%df_{(-%%d)}^{(+%%d)}', expo);
                    str = sprintf(form, num_disp, err_down_disp, err_up_disp);
                end
            else
                if err_up_disp == 0 && err_down_disp == 0
                    str = sprintf('%+*.0f', significant - expo, num_disp);
                elseif err_up_disp == err_down_disp
                    str = sprintf('%+*.0f(%d)', significant - expo, num_disp, err_up_disp*10^-expo);
                else
                    str = sprintf('%+*.0f_{(-%d)}^{(+%d)}', significant - expo, num_disp, err_down_disp*10^-expo, err_up_disp*10^-expo);
                end
            end
        end
    end
end

function [expo, err_disp, num_disp] = get_disp(significant, err, num)
    expo = log10(10^significant / err);
    if expo == floor(expo)
        expo = floor(expo) - 1;
    else
        expo = floor(expo);
    end
    err_disp = round(err*10^expo);
    num_disp = round(num*10^expo)*10^-expo;
end

function level = get_level(x)
    x = abs(x);
    if x == 0
        level = 0;
    elseif x < 1
        level = floor(log10(x));
    else
        level = floor(log10(x));
    end
end
