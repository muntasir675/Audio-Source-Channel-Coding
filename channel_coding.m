function channel_coding()
    % Load Stage 1 results
    if ~exist('stage1_results.mat','file')
        error('Stage2: run Stage 1 first to generate stage1_results.mat');
    end
    load('stage1_results.mat');

    % System parameters
    flip_probs_uncoded = [0 0.001 0.005 0.01 0.02];
    flip_probs_coded   = [0 0.01 0.02 0.03 0.05];

    constraint_length = 3;
    generators        = [7 5];
    trellis = poly2trellis(constraint_length, generators);
    tblen   = 5 * constraint_length;

    % Helper: symbol->audio
    symbol_to_audio = @(s) ((s - 33) ./ scale_factor);

    % Uncoded BSC (Huffman only)
    fprintf('\n--- Stage 2: Uncoded transmission (Huffman only) ---\n');
    for f = flip_probs_uncoded
        fprintf('   Testing uncoded flip probability f = %.3f ...\n', f);
        
        % Apply BSC errors
        Errors = rand(size(enc_bits)) < f;
        rx_bits = xor(enc_bits, Errors);
        
        % Huffman decode
        Nsym = numel(x_sym);
        dsym = zeros(Nsym,1);
        idx = 1;
        L = numel(rx_bits);
        incomplete = false;
        
        for n = 1:Nsym
            node = root;
            steps = 0;
            max_steps = 1000; % Prevent infinite loops
            
            while (nodes_left(node) ~= 0 || nodes_right(node) ~= 0) && steps < max_steps
                if idx > L
                    % Ran out of bits
                    dsym(n:end) = 33;
                    incomplete = true;
                    break;
                end
                
                b = rx_bits(idx);
                idx = idx + 1;
                steps = steps + 1;
                
                if b == 0
                    node = nodes_left(node);
                else
                    node = nodes_right(node);
                end
            end
            
            if incomplete
                break;
            end
            
            if steps >= max_steps
                % Corrupted data - fill rest with neutral
                dsym(n:end) = 33;
                fprintf('      WARNING: Decoding corrupted at symbol %d/%d\n', n, Nsym);
                break;
            end
            
            % Should be at leaf node
            if nodes_left(node) == 0 && nodes_right(node) == 0
                dsym(n) = nodes_sym(node);
            else
                dsym(n) = 33; % Shouldn't happen but safety
            end
        end
        
        x_uncoded = symbol_to_audio(dsym);
        fname = sprintf('stage2_uncoded_f_%0.3f.wav', f);
        audiowrite(fname, x_uncoded, Fs);
        fprintf('      -> Saved uncoded audio: %s\n', fname);
    end

    % Coded (convolutional) transmission
    fprintf('\n--- Stage 2: Coded transmission (conv + Huffman) ---\n');
    tx_bits = double(enc_bits(:));

    for f = flip_probs_coded
        fprintf('   Testing coded flip probability f = %.3f ...\n', f);

        % Convolutional encode (rate 1/2)
        coded_bits = convenc(tx_bits, trellis);
        coded_bits = coded_bits(:).';

        % BSC on coded bits
        Errors_coded = rand(size(coded_bits)) < f;
        rx_coded     = xor(coded_bits, Errors_coded);

        % Viterbi decode (hard decision)
        dec_conv = vitdec(double(rx_coded), trellis, tblen, 'trunc', 'hard');

        % Adjust to original Huffman length
        dec_conv = dec_conv(:).';
        if numel(dec_conv) < numel(enc_bits)
            dec_conv = [dec_conv, zeros(1, numel(enc_bits)-numel(dec_conv))];
        end
        dec_conv = dec_conv(1:numel(enc_bits));
        rx_bits2 = logical(dec_conv);

        % Huffman decode
        Nsym = numel(x_sym);
        dsym2 = zeros(Nsym,1);
        idx = 1;
        L2 = numel(rx_bits2);
        incomplete = false;
        
        for n = 1:Nsym
            node = root;
            steps = 0;
            max_steps = 1000;
            
            while (nodes_left(node) ~= 0 || nodes_right(node) ~= 0) && steps < max_steps
                if idx > L2
                    dsym2(n:end) = 33;
                    incomplete = true;
                    break;
                end
                
                b = rx_bits2(idx);
                idx = idx + 1;
                steps = steps + 1;
                
                if b == 0
                    node = nodes_left(node);
                else
                    node = nodes_right(node);
                end
            end
            
            if incomplete
                break;
            end
            
            if steps >= max_steps
                dsym2(n:end) = 33;
                fprintf('      WARNING: Decoding corrupted at symbol %d/%d\n', n, Nsym);
                break;
            end
            
            if nodes_left(node) == 0 && nodes_right(node) == 0
                dsym2(n) = nodes_sym(node);
            else
                dsym2(n) = 33;
            end
        end

        x_coded = symbol_to_audio(dsym2);
        fname2 = sprintf('stage2_coded_f_%0.3f.wav', f);
        audiowrite(fname2, x_coded, Fs);
        fprintf('      -> Saved coded audio: %s\n', fname2);
    end

    fprintf('\nStage 2 completed. Listen to the saved WAV files to judge acceptable flip probabilities.\n');
end
