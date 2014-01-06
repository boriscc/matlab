classdef map_2d < handle
% MAP_2D  Create 2D-surfs using few evaluations
%   MAP_2D generates data to use with e.g. SURF using as few evaluations of
%     FUN as possible, interpolating the other mesh points that are needed.
%   C = MAP_2D(PREFIX) loads a MAP_2D object from files with specified prefix
%   C = MAP_2D(PREFIX, FUN, X0, DX, X1, Y0, DY, Y1, ORDER = 3)
%     Initializes and stores a MAP_2D object.
%     FUN should be a function handle of type @(X, Y, PREFIX)
%     DX, DY should be as large as possible, but should not be larger than the
%       smallest feature of FUN that you want to resolve. The actual values
%       for DX and DY used by MAP_2D may be smaller than those specified.
%     ORDER is the number of points to initially use for each DX/DY length.
%       End points included. Must be an odd, positive integer >= 3.
%   C.calc(abstol, reltol, maxeval = Inf, maxlevel = Inf)
%     generate surface. Stop if tolerance is achieved or if maxeval is
%     reached. Do not allow more levels than maxlevel. The level is the
%     number of refinements of the grid that has been made. Initially it is 1.
%   C.interpolate() - get interpolated values for non-evaluted points
%     This is done automatically by calc.
%   C.set_max_interpolate_dxy(DX, DY)
    properties(GetAccess = 'public', SetAccess = 'private')
        prefix;
        f;
        order;
        s_tot;
        x0, x1, x_mul, x_tot;
        y0, y1, y_mul, y_tot;
        level_nr, level_reached;
        x, y, z;
        z_nr;
        z_calc;
        z_calc_nr;
        z_calc_time;
        err_max;
    end
    
    properties(GetAccess = 'public', SetAccess = 'public')
        save_after_calc;
    end
    
    methods(Access = private)
        function zi = interp1(~, z)
            zi = zeros(2 * size(z) - 1);
            zi(2:2:end) = 0.5 * (z(2:end) + z(1:(end-1)));
            zi(1:2:end) = z;
        end
        
        function zi = interp2(~, z)
            zi = zeros(2 * size(z) - 1);
            zi(1:2:end, 1:2:end) = z;
            zi(1:2:end, 2:2:end) = 0.5 * (zi(1:2:end, 1:2:(end-2)) + zi(1:2:end, 3:2:end));
            zi(2:2:end, 1:2:end) = 0.5 * (zi(1:2:(end-2), 1:2:end) + zi(3:2:end, 1:2:end));
            zi(2:2:end, 2:2:end) = 0.5 * (zi(2:2:end, 1:2:(end-2)) + zi(2:2:end, 3:2:end));
        end
        
        function precalc(obj)
            fprintf('--- Calculating first level ---\n');
            step = 2^(obj.level_nr-1);
            for ix = 1:step:obj.x_tot
                for iy = 1:step:obj.y_tot
                    obj.f_calc(ix, iy);
                end
            end
            fprintf('--- DONE Calculating first level ---\n');
        end
        
        function do_interpolation(obj, offx, offy, level_cur)
            if(level_cur == 0)
                return;
            end
            
            if(obj.x_tot == 1)
                n1x = 1;
            else
                n1x = obj.order;
            end
            if(obj.y_tot == 1)
                n1y = 1;
            else
                n1y = obj.order;
            end
            step1 = 2^level_cur;
            xidx1 = offx + 1 + step1*(0:(n1x-1));
            yidx1 = offy + 1 + step1*(0:(n1y-1));
            
            n2x = 2 * n1x - 1;
            n2y = 2 * n1y - 1;
            step2 = 2^(level_cur-1);
            xidx2 = offx + 1 + step2*(0:(n2x-1));
            yidx2 = offy + 1 + step2*(0:(n2y-1));
            
            if(~isempty(obj.z_calc(yidx2, xidx2) == 0))
                for iz = 1:obj.z_nr
                    if(numel(yidx1) == 1)
                        data = obj.interp1(obj.z(yidx1, xidx1, iz));
                    else
                        data = obj.interp2(obj.z(yidx1, xidx1, iz));
                    end
                    for nx = 1:length(xidx2);
                        ix = xidx2(nx);
                        for ny = 1:length(yidx2)
                            iy = yidx2(ny);
                            if(obj.z_calc(iy, ix) == 0)
                                obj.z(iy, ix, iz) = data(ny, nx);
                            end
                        end
                    end
                end
            end
            
            for ky = 1:2
                if(ky == 2 && n1y == 1)
                    break;
                end
                for kx = 1:2
                    obj.do_interpolation(offx + (kx - 1) * (n1x - 1) * 2^(level_cur-1), offy + (ky - 1) * (n1y - 1) * 2^(level_cur-1), level_cur - 1);
                end
            end
        end
        
        function z = f_calc(obj, ix, iy)
            if(obj.z_calc(iy, ix))
                z = obj.z(iy, ix);
            else
                t0 = tic;
                z = obj.f(obj.x(iy, ix), obj.y(iy, ix), obj.prefix);
                obj.z_calc_time = obj.z_calc_time + toc(t0);
                if(obj.z_calc_nr == 0 && numel(z) > 1)
                    obj.z = zeros(size(obj.z, 1), size(obj.z, 2), numel(z));
                    obj.z_nr = numel(z);
                end
                obj.z(iy, ix, :) = z;
                obj.z_calc(iy, ix) = 1;
                obj.z_calc_nr = obj.z_calc_nr + 1;
                if(obj.save_after_calc)
                    obj.save();
                end
            end
        end
        
        function err_max_ret = map_refine(obj, offx, offy, level_cur, abstol, reltol, tol_lim)
            if(obj.x_tot == 1)
                n1x = 1;
            else
                n1x = obj.order;
            end
            if(obj.y_tot == 1)
                n1y = 1;
            else
                n1y = obj.order;
            end
            step1 = 2^level_cur;
            xidx1 = offx + 1 + step1*(0:(n1x-1));
            yidx1 = offy + 1 + step1*(0:(n1y-1));

            n2x = 2 * n1x - 1;
            n2y = 2 * n1y - 1;
            step2 = 2^(level_cur-1);
            xidx2 = offx + 1 + step2*(0:(n2x-1));
            yidx2 = offy + 1 + step2*(0:(n2y-1));

            xidx1h = xidx1(1:2:end);
            yidx1h = yidx1(1:2:end);
            err_max_cur = -1;
            for iz = 1:obj.z_nr
                calc = obj.z(yidx1, xidx1, iz);
                if(numel(yidx1h) == 1)
                    inter = obj.interp1(obj.z(yidx1h, xidx1h, iz));
                else
                    inter = obj.interp2(obj.z(yidx1h, xidx1h, iz));
                end
                err_max_tmp = max(max(abs(calc - inter) ./ (abstol + reltol * abs(calc))));
                if(err_max_tmp > err_max_cur)
                    err_max_cur = err_max_tmp;
                end
            end
            
            if(err_max_cur < tol_lim)
                err_max_ret = err_max_cur;
                return;
            end

            if(obj.level_reached > level_cur)
                obj.level_reached = level_cur;
            end

            if(level_cur == 0)
                err_max_ret = err_max_cur;
                return;
            end

            for ix = xidx2
                for iy = yidx2
                    obj.f_calc(ix, iy);
                end
            end
            
            do_interp = ones(2, 2);
            err_max_sec = zeros(2, 2);
            
            for iz = 1:obj.z_nr
                calc = obj.z(yidx2, xidx2, iz);
                if(numel(yidx1) == 1)
                    inter = obj.interp1(obj.z(yidx1, xidx1, iz));
                else
                    inter = obj.interp2(obj.z(yidx1, xidx1, iz));
                end
                err = obj.get_error(abs(calc - inter) ./ (abstol + reltol * abs(calc)));

                for ky = 1:2
                    for kx = 1:2
                        if(err(ky, kx) > err_max_sec(ky, kx))
                            err_max_sec(ky, kx) = err(ky, kx);
                        end
                        if(err(ky, kx) >= tol_lim)
                            do_interp(ky, kx) = 0;
                        end
                    end
                end
            end
            err_max_ret = 0;
            for ky = 1:2
                if(ky == 2 && n1y == 1)
                    break;
                end
                for kx = 1:2
                    if(do_interp(ky, kx) == 1)
                        err_max_cur = err_max_sec(ky, kx);
                    else
                        err_max_cur = obj.map_refine(offx + (kx - 1) * (n1x - 1) * 2^(level_cur-1), offy + (ky - 1) * (n1y - 1) * 2^(level_cur-1), level_cur - 1, abstol, reltol, tol_lim);
                    end
                    if(err_max_cur > err_max_ret)
                        err_max_ret = err_max_cur;
                    end
                end
            end
        end

        function err = get_error(~, diff)
            n1 = (size(diff, 1) + 1) / 2;
            n2 = (size(diff, 2) + 1) / 2;
            err(1, 1) = max(max(diff(1:n1, 1:n2)));
            err(1, 2) = max(max(diff(1:n1, n2:end)));
            err(2, 1) = max(max(diff(n1:end, 1:n2)));
            err(2, 2) = max(max(diff(n1:end, n2:end)));
        end
        
        function change_size(obj, ix1, ix2, iy1, iy2)
            nx = ix2 - ix1 + 1;
            ny = iy2 - iy1 + 1;
            M = zeros(ny, nx);
            oix = max(1, ix1):min(size(obj.x, 2), ix2);
            oiy = max(1, iy1):min(size(obj.y, 1), iy2);
            nix = oix + 1 - ix1;
            niy = oiy + 1 - iy1;
            zz = zeros(ny, nx, obj.z_nr);
            for iz = 1:obj.z_nr
                zz(niy, nix, iz) = obj.z(oiy, oix, iz);
            end
            obj.z = zz;
            M(niy, nix) = obj.z_calc(oiy, oix);
            obj.z_calc = M;
            obj.z_calc_nr = sum(sum(obj.z_calc));
            [ obj.x, obj.y ] = meshgrid(linspace(obj.x0, obj.x1, obj.x_tot), linspace(obj.y0, obj.y1, obj.y_tot));
        end
    end
    
    methods
        function obj = map_2d(prefix, f, x0, dx, x1, y0, dy, y1, order)
            if(nargin == 1)
                obj.load(prefix);
            elseif(nargin >= 8)
                if(nargin < 9)
                    order = 3;
                end
                assert(order >= 3, 'Need order >= 3');
                assert(mod(order, 2) == 1, 'Need odd order');
                obj.prefix = prefix;
                obj.f = f;
                obj.order = order;
                if(x0 > x1)
                    tmp = x0;
                    x0 = x1;
                    x1 = tmp;
                end
                if(y0 > y1)
                    tmp = y0;
                    y0 = y1;
                    y1 = tmp;
                end
                obj.x0 = x0;
                obj.x1 = x1;
                obj.y0 = y0;
                obj.y1 = y1;
                if(x0 == x1)
                    obj.x_mul = 0;
                else
                    obj.x_mul = max(ceil(abs(x1 - x0) / ((obj.order-1) * dx)), 1);
                end
                if(y0 == y1)
                    obj.y_mul = 0;
                else
                    obj.y_mul = max(ceil(abs(y1 - y0) / ((obj.order-1) * dy)), 1);
                end
                obj.level_nr = 1;
                obj.level_reached = 1;
                obj.s_tot = 1 + 2^(obj.level_nr-1)*(obj.order-1);
                obj.x_tot = 1 + (obj.s_tot - 1) * obj.x_mul;
                obj.y_tot = 1 + (obj.s_tot - 1) * obj.y_mul;
                
                [ obj.x, obj.y ] = meshgrid(linspace(obj.x0, obj.x1, obj.x_tot), linspace(obj.y0, obj.y1, obj.y_tot));
                obj.z = zeros(obj.y_tot, obj.x_tot);
                obj.z_nr = 1;
                obj.z_calc = zeros(obj.y_tot, obj.x_tot);
                obj.z_calc_nr = 0;
                obj.err_max = Inf;
                
                obj.save_after_calc = 0;
                obj.save();
            else
                error('Bad number of arguments');
            end
        end
    end
    
    methods(Access = public)
        function load(obj, prefix)
            obj.prefix = prefix;
            filename = sprintf('%s.func', obj.prefix);
            fp = fopen(filename, 'r');
            assert(fp ~= -1, sprintf('Could not open func-file "%s"', filename));
            obj.f = str2func(fscanf(fp, '%s'));
            fclose(fp);
            
            data = num2cell(dlmread(sprintf('%s.param', obj.prefix)));
            [ obj.order, obj.s_tot, obj.x0, obj.x1, obj.x_mul, obj.x_tot, obj.y0, obj.y1, obj.y_mul, obj.y_tot, obj.level_nr, obj.level_reached, obj.z_nr, obj.z_calc_nr, obj.err_max, obj.save_after_calc ] = deal(data{:});
            obj.x = dlmread(sprintf('%s.x', obj.prefix));
            obj.y = dlmread(sprintf('%s.y', obj.prefix));
            obj.z = zeros(size(obj.x, 1), size(obj.x, 2), obj.z_nr);
            for n = 1:obj.z_nr
                obj.z(:, :, n) = dlmread(sprintf('%s.z_%d', obj.prefix, n));
            end
            obj.z_calc = dlmread(sprintf('%s.z_calc', obj.prefix));
        end
        
        function save(obj, new_prefix)
            if(nargin > 1)
                obj.prefix = new_prefix;
            end
            filename = sprintf('%s.func', obj.prefix);
            fp = fopen(filename, 'w');
            assert(fp ~= -1, sprintf('SAVE FAILED: Could not open func-file "%s" for writing.', filename));
            fprintf(fp, '%s', func2str(obj.f));
            fclose(fp);
            
            dlmwrite(sprintf('%s.param', obj.prefix), [ obj.order, obj.s_tot, obj.x0, obj.x1, obj.x_mul, obj.x_tot, obj.y0, obj.y1, obj.y_mul, obj.y_tot, obj.level_nr, obj.level_reached, obj.z_nr, obj.z_calc_nr, obj.err_max, obj.save_after_calc ]);
            dlmwrite(sprintf('%s.x', obj.prefix), obj.x);
            dlmwrite(sprintf('%s.y', obj.prefix), obj.y);
            for n = 1:obj.z_nr
                dlmwrite(sprintf('%s.z_%d', obj.prefix, n), obj.z(:, :, n));
            end
            dlmwrite(sprintf('%s.z_calc', obj.prefix), obj.z_calc);
            filename = sprintf('%s.txt', obj.prefix);
            fp = fopen(filename, 'w');
            assert(fp ~= -1, sprintf('Could not open "%s" for writing, no .txt-file saved (but everything else is saved).', filename));
            for ix = 1:obj.x_tot
                for iy = 1:obj.y_tot
                    fprintf(fp, '%g %g', obj.x(iy, ix), obj.y(iy, ix));
                    for iz = 1:obj.z_nr
                        fprintf(fp, ' %g', obj.z(iy, ix, iz));
                    end
                    fprintf(fp, '\n');
                end
            end
            fclose(fp);
        end
        
        function set_xy01(obj, x0, x1, y0, y1)
            if(x0 > x1)
                tmp = x0;
                x0 = x1;
                x1 = tmp;
            end
            if(y0 > y1)
                tmp = y0;
                y0 = y1;
                y1 = tmp;
            end
            
            if(x0 == x1)
                assert(obj.x_tot == 1);
                dx = 0;
                n1x = 0;
                n2x = 0;
            else
                dx = obj.x(1, obj.s_tot) - obj.x(1, 1);
                n1x = ceil((obj.x0 - x0) / dx);
                n2x = ceil((x1 - obj.x1) / dx);
                obj.x_mul = obj.x_mul + n1x + n2x;
            end
            if(y0 == y1)
                assert(obj.y_tot == 1);
                dy = 0;
                n1y = 0;
                n2y = 0;
            else
                dy = obj.y(obj.s_tot, 1) - obj.y(1, 1);
                n1y = ceil((obj.y0 - y0) / dy);
                n2y = ceil((y1 - obj.y1) / dy);
                obj.y_mul = obj.y_mul + n1y + n2y;
            end
            
            x_tot_old = obj.x_tot;
            y_tot_old = obj.y_tot;
            obj.x_tot = 1 + (obj.s_tot - 1) * obj.x_mul;
            obj.y_tot = 1 + (obj.s_tot - 1) * obj.y_mul;
            m = obj.s_tot - 1;
            obj.x0 = obj.x0 - n1x * dx;
            obj.x1 = obj.x1 + n2x * dx;
            obj.y0 = obj.y0 - n1y * dy;
            obj.y1 = obj.y1 + n2y * dy;
            obj.change_size(1 - n1x*m, x_tot_old + n2x*m, 1 - n1y*m, y_tot_old + n2y*m);
        end
        
        function set_dx(obj, dx)
            x_mul_needed = max(ceil(abs(obj.x1 - obj.x0) / ((obj.order-1) * dx)), 1);
            
            while(x_mul_needed > obj.x_mul)
                obj.x_mul = obj.x_mul * 2;
                obj.x_tot = 1 + (obj.s_tot - 1) * obj.x_mul;
                [ obj.x, obj.y ] = meshgrid(linspace(obj.x0, obj.x1, obj.x_tot), linspace(obj.y0, obj.y1, obj.y_tot));
                M = zeros(size(obj.x));
                zz = zeros(size(obj.x, 1), size(obj.x, 2), obj.z_nr);
                for iz = 1:obj.z_nr
                    zz(:, 1:2:end, iz) = obj.z(:, :, iz);
                end
                obj.z = zz;
                M(:, 1:2:end) = obj.z_calc;
                obj.z_calc = M;
                obj.z_calc_nr = sum(sum(obj.z_calc));
            end
        end
        
        function set_dy(obj, dy)
            y_mul_needed = max(ceil(abs(obj.y1 - obj.y0) / ((obj.order-1) * dy)), 1);
            
            while(y_mul_needed > obj.y_mul)
                obj.y_mul = obj.y_mul * 2;
                obj.y_tot = 1 + (obj.s_tot - 1) * obj.y_mul;
                [ obj.x, obj.y ] = meshgrid(linspace(obj.x0, obj.x1, obj.x_tot), linspace(obj.y0, obj.y1, obj.y_tot));
                M = zeros(size(obj.x));
                zz = zeros(size(obj.x, 1), size(obj.x, 2), obj.z_nr);
                for iz = 1:obj.z_nr
                    zz(1:2:end, :, iz) = obj.z(:, :, iz);
                end
                obj.z = zz;
                M(1:2:end, :) = obj.z_calc;
                obj.z_calc = M;
                obj.z_calc_nr = sum(sum(obj.z_calc));
            end
        end
        
        function add_level(obj)
            obj.s_tot = 2 * obj.s_tot - 1;
            obj.x_tot = 2 * obj.x_tot - 1;
            obj.y_tot = 2 * obj.y_tot - 1;
            obj.level_nr = obj.level_nr + 1;
            obj.level_reached = obj.level_reached + 1;
            [ obj.x, obj.y ] = meshgrid(linspace(obj.x0, obj.x1, obj.x_tot), linspace(obj.y0, obj.y1, obj.y_tot));
            M = zeros(size(obj.x));
            zz = zeros(size(obj.x, 1), size(obj.x, 2), obj.z_nr);
            for iz = 1:obj.z_nr
                zz(1:2:end, 1:2:end, iz) = obj.z(:, :, iz);
            end
            obj.z = zz;
            M(1:2:end, 1:2:end) = obj.z_calc;
            obj.z_calc = M;
            obj.z_calc_nr = sum(sum(obj.z_calc));
        end
        
        function set_max_interpolate_dxy(obj, dx, dy)
            if(obj.y_mul == 0)
                assert(dx > 0, 'Need dx > 0');
                while(abs(obj.x(1, 1) - obj.x(1, 2)) > dx)
                    obj.add_level();
                end
            elseif(obj.x_mul == 0)
                assert(dy> 0, 'Need dy > 0');
                while(abs(obj.y(1, 1) - obj.y(2, 1)) > dy)
                    obj.add_level();
                end
            else
                assert(dx > 0, 'Need dx > 0');
                assert(dy> 0, 'Need dy > 0');
                while(abs(obj.x(1, 1) - obj.x(1, 2)) > dx && abs(obj.y(1, 1) - obj.y(2, 1)) > dy)
                    obj.add_level();
                end
            end
        end
        
        function interpolate(obj)
            fprintf('Interpolating... ');
            if(obj.y_mul == 0)
                ic = find(obj.z_calc == 1);
                ii = find(obj.z_calc == 0);
                for n = 1:obj.z_nr
                    obj.z(ii + (n-1) * numel(obj.x)) = interp1(obj.x(ic), obj.z(ic + (n-1) * numel(obj.x)), obj.x(ii));
                end
            else
                for n1 = 0:max(obj.x_mul - 1, 0)
                    for n2 = 0:max(obj.y_mul - 1, 0)
                        %fprintf('Interpolating block (%d, %d)\n', n1, n2);
                        obj.do_interpolation(n1 * (obj.s_tot - 1), n2 * (obj.s_tot - 1), obj.level_nr - 1);
                    end
                end
            end
            fprintf('done.\n');
        end
        
        function calc(obj, abstol, reltol, eval_max, level_max)
            if(nargin < 5)
                level_max = Inf;
            end
            if(nargin < 4)
                eval_max = Inf;
            end
            assert(abstol >= 0, 'Need abstol >= 0');
            assert(reltol >= 0, 'Need reltol >= 0');
            assert(abstol > 0 || reltol > 0, 'Need one of abstol and reltol > 0');
            z_calc_pre = obj.z_calc_nr;
            obj.z_calc_time = 0;
            t0 = tic;
            obj.precalc();

            tol_lim = Inf;
            obj.err_max = Inf;
            while obj.err_max > 1 && obj.z_calc_nr < z_calc_pre + eval_max
                obj.err_max = -1;
                for n1 = 0:max(obj.x_mul - 1, 0)
                    for n2 = 0:max(obj.y_mul - 1, 0)
                        err_max_cur = obj.map_refine(n1 * (obj.s_tot - 1), n2 * (obj.s_tot - 1), obj.level_nr - 1, abstol, reltol, tol_lim);
                        obj.err_max = max(obj.err_max, err_max_cur);
                    end
                end
                fprintf('eval_nr = %d/%d, err_max = %.3e\n', obj.z_calc_nr, numel(obj.x), obj.err_max);
                tol_lim = obj.err_max;
                if(obj.level_reached == 0 || obj.z_calc_nr == numel(obj.x))
                    if(obj.level_nr < level_max)
                        fprintf('Adding a level... ');
                        obj.add_level();
                        fprintf('done\n');
                        obj.precalc();
                        tol_lim = Inf;
                        obj.err_max = Inf;
                    end
                end
            end
            obj.interpolate();
            obj.save();
            time_tot = toc(t0);
            fprintf('Time spent evaluating f: %.1f%%\n', 100 * obj.z_calc_time / time_tot);
        end
    end
end

