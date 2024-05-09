simulation()

function simulation()
    symbol_duration = 1; 
    number_of_bits = 1500000;
    transmitted_bits = randi([0,1],1, number_of_bits);
    M = [2,8];
    SNR = linspace(0,30,31);
    ser_arr = zeros(2,31); % symbol error rate array
    theoretical_Pe = zeros(2,31); % analytical Pe values

    for i = 1:2
        simulated_M = M(i);
        
        for j = 1:31

            simulated_SNR = SNR(j);
            
            theoretical_Pe(i,j) = analytic_err_calculator(simulated_M, simulated_SNR);
        
            [modulated_symbol,number_of_symbols] = m_ask_modulator(transmitted_bits, simulated_M, 5, symbol_duration);

            [received_bits] = m_ask_demodulator(modulated_symbol, simulated_M, 5, symbol_duration,number_of_symbols, simulated_SNR);
            
            ser_arr(i,j) = ser_calculator(transmitted_bits, received_bits, simulated_M, number_of_symbols);
        
        end

    end
    
    semilogy(SNR, theoretical_Pe(1,:), '-x',SNR,theoretical_Pe(2,:), '-x', SNR, ser_arr(1,:), '-o', SNR, ser_arr(2,:), '-o');
    ylim([10^-6, 10^0])
    xlim([0,30])
    ylabel('Symbol Error Rate (SER)')
    xlabel('SNR (dB)')
    legend('M = 2 Theoretical', 'M = 8 Theoretical', 'M = 2 Simulated', 'M = 8 Simulated')
    
end

function [Pe] = analytic_err_calculator(M, SNR)

    coeffs = 0:3:3*(M-1);
    
    E_s = mean(coeffs.*coeffs) / 2; % Average symbol energy

    linear_SNR = 10 .^(SNR / 10);

    N0 = E_s ./ linear_SNR;
    
    Pe = 2 * (1 - 1 / M) * qfunc(sqrt(E_s * 3 / ((M-1) * (2*M-1) * N0)));
end

function [ser] = ser_calculator(transmitted_bits, received_bits, M, number_of_symbols)

    bits_per_symbol = log2(M);
    correct_symbol_count = 0;
    false_symbol_count = 0;

    for i = 1:number_of_symbols

        if transmitted_bits(bits_per_symbol * (i-1) + 1 : bits_per_symbol * i) == received_bits(bits_per_symbol * (i-1) + 1 : bits_per_symbol * i)
            correct_symbol_count = correct_symbol_count + 1;
        else
            false_symbol_count = false_symbol_count + 1;
        end

    end

    ser = (false_symbol_count / number_of_symbols);
end

function [received_bits] = m_ask_demodulator(received_signal, M, frequency, symbol_duration, number_of_symbols, snr)

    symbol_time_domain = linspace(0,symbol_duration);

    symbol_duration_in_vector = length(received_signal) / number_of_symbols; % the number of elements in the vector corresponding to one symbol

    ref_symbol_energy = 9/2; % the energy of the symbol in the correlated receiver

    bits_per_symbol = log2(M);
    
    s_1_t = 3 * cos(2 * pi * frequency * symbol_time_domain); % the symbol in the correlated receiver

    z_k_T = zeros(1,number_of_symbols); % the values getting fed to the decision circuit

    received_bits = zeros(1, bits_per_symbol * number_of_symbols);

    E_s= (M-1) * (2*M - 1) * ref_symbol_energy^2 / 12; % E_s

    E_s_dB = 10 * log10(E_s); % E_s in dB

    index = 1;

    for i = 1: number_of_symbols
        corresponding_symbol = received_signal(symbol_duration_in_vector * (i - 1) + 1 : symbol_duration_in_vector * i);
        
        multiplied_function = corresponding_symbol .* s_1_t; % the output of the multiplication in correlated receiver

        z_k_T(index) = multiplied_function(1) / 2; % getting the related value from the array instead of integrating 

        index = index + 1;
    end
    
    z_k_T = awgn(z_k_T, snr, E_s_dB, 'dB'); % Adding noise

    % Decision Circuit Algorithm
    for i = 1:number_of_symbols
        corresponding_bit_sequence =  zeros(1,bits_per_symbol); 
        corresponding_voltage = z_k_T(i);
        if M == 8
            if corresponding_voltage < 2.25
                corresponding_bit_sequence = [0 0 0];
            
            elseif (2.25 <= corresponding_voltage) && (corresponding_voltage < 6.75)
                corresponding_bit_sequence = [0 0 1];
              
            elseif (6.75 <= corresponding_voltage) && (corresponding_voltage < 11.25)
                corresponding_bit_sequence = [0 1 0];
    
            elseif (11.25 <= corresponding_voltage) && (corresponding_voltage < 15.75)
                corresponding_bit_sequence = [0 1 1];
            
            elseif (15.75 <= corresponding_voltage) && (corresponding_voltage < 20.25)
                corresponding_bit_sequence = [1 0 0];
    
            elseif (20.25 <= corresponding_voltage) && (corresponding_voltage < 24.75)
                corresponding_bit_sequence = [1 0 1];
    
            elseif (24.75 <= corresponding_voltage) &&( corresponding_voltage < 29.25)
                corresponding_bit_sequence = [1 1 0];
        
            elseif 29.25 <= corresponding_voltage 
                corresponding_bit_sequence = [1 1 1];
            end
            received_bits(bits_per_symbol * (i-1) + 1 : bits_per_symbol * i) = corresponding_bit_sequence;
        
        elseif M == 2
            if corresponding_voltage < 2.25
                corresponding_bit_sequence = 0;
            
            elseif (2.25 <= corresponding_voltage) && (corresponding_voltage < 6.75)
                corresponding_bit_sequence = 1;
            
            end
            received_bits(bits_per_symbol * (i-1) + 1 : bits_per_symbol * i) = corresponding_bit_sequence;
            
        end

    end
end

function [modulated_symbol,number_of_symbols] = m_ask_modulator(bit_array, M, frequency, symbol_duration)

    array_length = length(bit_array);
    
    bits_per_symbol = log2(M);

    number_of_symbols = array_length / log2(M);

    final_symbol_length = 100 * symbol_duration;
    
    modulated_symbol = zeros(1,final_symbol_length);
    
    for i = 1:number_of_symbols
        
        symbol_time_domain = linspace((i-1) * symbol_duration, i * symbol_duration, 100);

        bit_sequence = bit_array( (i-1)*bits_per_symbol + 1 : i*bits_per_symbol);

        corresponding_decimal = binary_to_decimal(bit_sequence);

        amplitude = corresponding_decimal * 3; % in Volts

        symbol = amplitude * cos(2 * pi * frequency * symbol_time_domain);

        modulated_symbol((i-1) * 100 + 1: i * 100) =  symbol;
       
        
    end
end

function [decimal] =  binary_to_decimal(array)

    array_length = length(array);

    power = zeros(1,array_length - 1);

    for i = 0 : array_length - 1
        power(i + 1) = 2^((array_length - i) - 1);
    end
    
    decimal = power * transpose(array);

end