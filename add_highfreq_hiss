function [x_noisy, hiss] = add_hiss(x, snr_dB)


    if size(x,2) > 1
        x = mean(x, 2); 
    end

    rng(0);

    sigPower = mean(x.^2);

    noisePower = sigPower / (10^(snr_dB/10));

    hiss = sqrt(noisePower) * randn(size(x));

    x_noisy = x + hiss;
end
