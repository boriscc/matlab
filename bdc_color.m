function C = bdc_color(n, c1, c2, power, non_saturated, mode)
% C = BDC_COLOR(N, C1, C2, POWER=1.0, NON_SATURATED=1.0, MODE='linear') Generate color-matrices for colormap
%   N: Number of colors
%   C1: Starting color, 3-vector, elements in range [0, 1]
%   C2: Ending color, 3-vector, elements in range [0, 1]
%   POWER (> 0): != 1 => C -> c1 + (c2 - c1) * ((C - c1) / (c2 - c1)).^power
%   NON_SATURATED: Must be in [-1, 1]. If > 0, then that much of the lower part
%                  will be non saturated, the upper part will be saturated.
%                  If < 0, then that much of the upper part will be non saturated,
%                  the lower part will be saturated.
%   MODE: Must be 'linear'
    if(nargin < 6)
        mode = 'linear';
    end
    if(nargin < 5)
        non_saturated = 1;
    end
    if(nargin < 4)
        power = 1.0;
    end
    assert(numel(power) == 1, 'power must be a scalar');
    assert(power > 0, 'Must have power > 0');
    assert(numel(c1) == 3, 'c1 should be a 3-vector');
    assert(numel(c2) == 3, 'c2 should be a 3-vector');
    assert(numel(n) == 1, 'n must be a scalar');
    assert(n > 0, 'Must have n > 0');
    assert(numel(non_saturated) == 1, 'non_saturated must be a scalar');
    assert(abs(non_saturated) <= 1, 'Must have abs(non_saturated) <= 1');
    assert(min(c1) >= 0 && max(c1) <= 1, 'Elements of c1 must be in range [0, 1]');
    assert(min(c2) >= 0 && max(c2) <= 1, 'Elements of c2 must be in range [0, 1]');
    C = zeros(n, 3);
    n1 = abs(round(n * non_saturated));
    for m = 1:3
        if(c2(m) == c1(m))
            C(:, m) = c1(m);
        else
            if non_saturated > 0
                if mode == 'linear'
                    C(1:n1, m) = linspace(c1(m), c2(m), n1);
                else
                    error('Invalid mode');
                end
                C((1+n1):end, m) = c2(m);
            else
                C(1:(end - n1), m) = c1(m);
                if mode == 'linear'
                    C((end - n1 + 1):end, m) = linspace(c1(m), c2(m), n1);
                else
                    error('Invalid mode');
                end
            end
        end
    end

    if(power ~= 1.0)
        C = c1 + (c2 - c1) .* ((C - c1) ./ (c2 - c1)).^power;
    end
end

