function [x_noisy, hum] = add_hum(x, Fs, f_hum, hum_amp_scale)

    if size(x,2) > 1
        x = mean(x, 2);
    end

    if nargin < 3
        f_hum = 60;   
    end
    if nargin < 4
        hum_amp_scale = 0.1; 
    end

    N = length(x);
    t = (0:N-1)'/Fs;

    A_hum = hum_amp_scale * max(abs(x));  
    hum   = A_hum * sin(2*pi*f_hum*t);

    x_noisy = x + hum;
end
