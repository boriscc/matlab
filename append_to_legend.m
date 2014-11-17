function append_to_legend(h, name, varargin)
% APPEND_TO_LEGEND  Append a new entry to a legend
%   APPEND_TO_LEGEND(plot_handle, legend_entry_name)
%   If there is no legend, a legend will be created.
%
% See also LEGEND PLOT
    [~, ~, outh, outm] = legend;
    if(isempty(outm))
        legend(h, name, varargin{:});
    else
        legend([outh; h], {outm{:}, name}, varargin{:});
    end
end

