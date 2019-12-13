close all;
clear all;
rand('seed', sum(100*clock));


%% variables
hist_res = 1000;
samples = 1000000;

%% define r.v. Y
x_sample = linspace(0,log(4),samples);
x1 = exp(-x_sample/2);
x2 = exp(-x_sample/2);
y = log(randsample(x1,samples)) - log(randsample(x2,samples));
y = y.*2;

%% plot the histogram
h_y = histogram(y,hist_res);
title('HISTOGRAM OF Y');
y_power = calculate_pow(y)



%% take the uniform quantization of it
y_uni_quantized = quantize_uniform(y, log(4), -log(4), 32);

%% calculate the uniform quantization error
uni_q_e = y - y_uni_quantized;
h_uni_q_e = histogram(uni_q_e, hist_res);
title('HISTOGRAM OF UNIFORM QUANTIZATION ERROR');
uni_mse = calculate_pow(uni_q_e)

%% SQNR and D in db
sqnr_db = 10*log10(y_power/uni_mse)
D_db = 10*log10(uni_mse)

%% Lloyd-Max Algorithm
boundaries = [linspace(-log(4),-0.4,10), linspace(-0.35,0.35,13), linspace(0.4,log(4),10)];
reconstruction = zeros(1,32);
mse_prev = 1;
mse_new = 1;

while (1)
    [hist_values, edges] = histcounts(y,1000);
    %define the reconstruction levels
    for x = 1:32
        total_data_num = 0;
        expectation_t = 0;
        i = 1;
        while (true)
            temp_num = (edges(i)+edges(i+1))/2;
            if ( temp_num >= boundaries(x) && temp_num < boundaries(x+1)  )
                total_data_num = total_data_num + hist_values(i);
                expectation_t = expectation_t + hist_values(i)*temp_num;
                if i < 1000
                    i = i + 1;
                else
                    break;
                end;
            else
                if total_data_num == 0 && i < 1000
                    i = i + 1;
                    continue;
                else
                    break;
                end;
            end;
        end;
        reconstruction(x) = expectation_t/total_data_num;
    end;
    %define new boundaries
    for x = 1:31
        boundaries(x+1) = (1/2)*(reconstruction(x) + reconstruction(x+1));
    end;
    %check new error 
    y_nonuni_quantized = quantize_nonuniform(y, boundaries, reconstruction);
    err_nonuni = y - y_nonuni_quantized;
    mse_new = calculate_pow(err_nonuni);
    %finishing condition (no more than %0.01 improvement on err)
    if( abs((mse_new-mse_prev)/mse_new) <= 0.0001)
        break;
    end;
    mse_prev = mse_new;
end;

mse_new;
sqnr_nonuni_db = 10*log10(y_power/mse_new)
improvement = abs((mse_new-uni_mse)/mse_new*100)
boundaries
reconstruction
histogram(err_nonuni,1000);

temp_record_data =  audiorecorder(200000,16,1);
disp('record')
recordblocking(temp_record_data, 5);
disp('record over')

recoded_sound = getaudiodata(temp_record_data, 'double');
recoded_sound = transpose(recoded_sound);
voicehist = histogram(recoded_sound,1000,'Normalization','pdf');

voicepower = calculate_pow(recoded_sound);
voicepower_in_dB = 10*log10(voicepower)

voice_uni_quantized = quantize_uniform(recoded_sound, 1, -1, 64);
err = recoded_sound - voice_uni_quantized;
histogram(err,1000,'Normalization','pdf') %pdf of error
err_power = calculate_pow(err);
sqnr = voicepower/err_power;
sqnr_db = 10*log10(sqnr)

%% mu compand
compressed = log(1+64.*abs(recoded_sound))/log(1+64).*sign(recoded_sound);
voice_nonuni_quantized = quantize_uniform(compressed, max(compressed), min(compressed), 64);
expanded = sign(compressed).*1/64.*((1+64).^abs(compressed)-1);
err_64 = recoded_sound - expanded;
err_power_64 = calculate_pow(err_64);
sqnr_64 = voicepower/err_power_64;
sqnr_db_64 = 10*log10(sqnr_64)



%% function for calculating power
function powerr = calculate_pow(fnc_in)
    [hist_values_e, edges] = histcounts(fnc_in,1000,'Normalization','pdf');
    total_e = 0;
    i = 0;
    x = linspace(min(edges), max(edges), 1000);
    width = x(2) - x(1);
    for k = x
        i = i + 1;
        total_e = total_e + (k^2)*width*hist_values_e(i);
    end;
    powerr = total_e;
end

%% function for non-uniform quantization
function y_quantized = quantize_nonuniform(y, boundaries, reconstruction)
    y_quantized = zeros(1,1000000);
    for x = 1:1000000
        for i = 1:32
            if (y(x) >= boundaries(i) && y(x) < boundaries(i+1) )
                y_quantized(x) = reconstruction(i);
                break;
            end
        end
    end
    return;
end

%% function for uniform quantization
function y_uni_quantized = quantize_uniform(y, y_max, y_min, level)
    uni_q_delt = (y_max-y_min)/level;
    y_uni_quantized = zeros(1, 1000000);
    for x = 1:1000000
        which_width = fix((y(x)+y_max)/uni_q_delt);
        quantized_value = y_min + (uni_q_delt/2) + which_width*uni_q_delt;
        y_uni_quantized(x) = quantized_value;
    end
    return;
end
