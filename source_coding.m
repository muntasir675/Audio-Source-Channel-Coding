function source_coding()
    % Parameters
    bits_per_sample = 6;                 % K
    audio_file      = 'sample.wav';

    % 1) Read and quantize audio
    if ~isfile(audio_file)
        error('Audio file "%s" not found in current folder.', audio_file);
    end
    [raw_signal, Fs] = audioread(audio_file);
    raw_signal = raw_signal(:,1);               % mono
    scale_factor = 2^(bits_per_sample-1);      % leave 1 bit for sign
    q_int       = floor(raw_signal * scale_factor); % approx [-32..31]
    x_sym       = q_int + 33;                  % map to [1..64]
    x_sym       = min(max(x_sym,1), 2^bits_per_sample); % safety clamp

    M = 2^bits_per_sample;                     % 64
    N = numel(x_sym);

    % 2) Original fixed-length bitstream using dec2bin
    bin_str        = dec2bin(x_sym, bits_per_sample);   % N x 6
    orig_bits_char = reshape(bin_str.', 1, []);         % serialize char
    orig_bits      = (orig_bits_char == '1');          % logical vector
    L_orig         = numel(orig_bits);

    % 3) PMF and entropy
    counts = histcounts(x_sym, 0.5:1:(M+0.5));
    pmf    = counts / N;
    figure; stem(1:M, pmf, 'filled'); xlabel('Symbol index'); ylabel('Probability');
    title('PMF of quantized audio symbols'); grid on;

    p_nonzero = pmf(pmf>0);
    H         = -sum(p_nonzero .* log2(p_nonzero));     % entropy
    N_min_bits = H * N;

    % 4) Custom Huffman dictionary (tree build) using arrays
    symbols = (1:M).';
    valid_idx = find(pmf>0);
    valid_sym = symbols(valid_idx);
    valid_p   = pmf(valid_idx);
    K_valid   = numel(valid_sym);
    

    % initialize nodes arrays
    nodes_p     = valid_p(:);                  % probabilities
    nodes_left  = zeros(K_valid,1);            % 0 means no child
    nodes_right = zeros(K_valid,1);
    nodes_sym   = valid_sym(:);                % leaf symbol or 0 for internal

    active = (1:K_valid).';
    % grow arrays as we add parents
    while numel(active) > 1
        % sort active indices by probability
        [~, order] = sort(nodes_p(active));
        active = active(order);

        a = active(1);
        b = active(2);

        new_idx = numel(nodes_p) + 1;
        nodes_p(new_idx,1) = nodes_p(a) + nodes_p(b);
        nodes_left(new_idx,1)  = a;
        nodes_right(new_idx,1) = b;
        nodes_sym(new_idx,1)   = 0; % internal

        % remove first two active, add parent
        active = [new_idx; active(3:end)];
    end

    root = active; % index of root node

    % 5) Get codes by iterative DFS (non-recursive) producing numeric vectors
    codes_valid = cell(K_valid,1);             % codes for valid_sym order
    % stack holds pairs {node_index, code_array}
    stack = {root, zeros(0,1)};                % begin with root and empty code (as column)

    while ~isempty(stack)
        item = stack(end,:);
        stack(end,:) = [];                     % pop
        node = item{1};
        prefix = item{2};

        if nodes_left(node) == 0 && nodes_right(node) == 0
            s = nodes_sym(node);
            k = find(valid_sym == s, 1);
            codes_valid{k} = prefix(:).';     % store as row vector
        else
            % push right then left so left is processed first (for lexicographic)
            if nodes_right(node) ~= 0
                stack(end+1,:) = {nodes_right(node), [prefix; 1]}; %#ok<AGROW>
            end
            if nodes_left(node) ~= 0
                stack(end+1,:) = {nodes_left(node), [prefix; 0]}; %#ok<AGROW>
            end
        end
    end

    % 6) Full dictionary 1..M
    dict_custom = cell(M,2);
    for i = 1:M
        dict_custom{i,1} = i;
        dict_custom{i,2} = zeros(0,1); % empty code for zero-prob symbols
    end
    for k = 1:K_valid
        s = valid_sym(k);
        dict_custom{s,2} = codes_valid{k}(:).'; % row vector 0/1
    end

    % 7) Encode using dictionary (efficient preallocation)
    % compute length per symbol for all samples
    lens_per_sym = zeros(N,1);
    for i=1:N
        lens_per_sym(i) = numel(dict_custom{x_sym(i),2});
    end
    totalLen = sum(lens_per_sym);
    enc_bits = false(1, totalLen);  % preallocate logical array
    ptr = 1;
    for n = 1:N
        code = dict_custom{x_sym(n),2};        % numeric row vector of 0/1
        Lc = numel(code);
        if Lc>0
            enc_bits(ptr:ptr+Lc-1) = logical(code);
            ptr = ptr + Lc;
        end
    end
    L_huff = numel(enc_bits);
    CR_huff = L_orig / L_huff;
    CR_theory = bits_per_sample / H;

    % 8) Custom decoder using tree (walk bits)
    dec_sym = zeros(N,1);
    idx = 1;
    L_enc = L_huff;
    for n = 1:N
        node = root;
        while nodes_left(node) ~= 0 || nodes_right(node) ~= 0
            if idx > L_enc
                error('Stage1: ran out of bits while decoding.');
            end
            b = enc_bits(idx);
            idx = idx + 1;
            if b == 0
                node = nodes_left(node);
            else
                node = nodes_right(node);
            end
        end
        dec_sym(n) = nodes_sym(node);
    end

    if ~isequal(dec_sym(:), x_sym(:))
        error('Stage1: decoded symbols do not match original.');
    end
% Example: print first 10 non-empty codes to Command Window
fprintf('\nSample of Huffman dictionary (symbol : code):\n');
count = 0;
for i = 1:M
    if ~isempty(dict_custom{i,2})
        code_str = num2str(dict_custom{i,2});
        code_str(code_str==' ') = [];   % remove spaces
        fprintf('%2d : %s\n', i, code_str);
        count = count + 1;
        if count == 10
            break;
        end
    end
end


    % 9) Reconstruct audio and save
    q_int_rec   = dec_sym - 33;
    x_rec       = q_int_rec / scale_factor;
    audiowrite('stage1_reconstructed_custom.wav', x_rec, Fs);



    % 10) Bonus: compare with built-in Huffman (if Communications Toolbox present)
    has_builtin = exist('huffmandict','file')==2 && exist('huffmanenco','file')==2;

    if has_builtin && N < 100000  % Add size check
        sym_cells = num2cell((1:M));
        dict_builtin = huffmandict(sym_cells, pmf);
        x_cells = num2cell(x_sym(:).');
        bits_builtin = huffmanenco(x_cells, dict_builtin);
        L_builtin = numel(bits_builtin);
        CR_builtin = L_orig / L_builtin;
        
        dec_cells = huffmandeco(bits_builtin, dict_builtin);
        dec_sym_bi = cell2mat(dec_cells);
        if ~isequal(dec_sym_bi(:), x_sym(:))
            warning('Stage1: built-in Huffman does not perfectly reconstruct.');
        end
    elseif has_builtin
        fprintf(' (Built-in comparison skipped: file too large)\n');
        L_builtin = NaN;
        CR_builtin = NaN;
    else
        L_builtin = NaN;
        CR_builtin = NaN;
    end


% Choose a later starting point so you avoid the initial silence
start_sample = Fs * 5;             % start after 5 second (for example)
Nplot        = min(5000, numel(raw_signal) - start_sample);
idx          = start_sample + (0:Nplot-1);
t            = (0:Nplot-1) / Fs;

figure;
subplot(2,1,1);
plot(t, raw_signal(idx));
xlabel('Time (s)');
ylabel('Amplitude');
title('Original audio waveform (segment)');
grid on;

subplot(2,1,2);
plot(t, x_rec(idx));
xlabel('Time (s)');
ylabel('Amplitude');
title('Huffman-decoded audio waveform (segment)');
grid on;


    % 11) Save variables for Stage 2
    save('stage1_results.mat', 'enc_bits', 'x_sym', 'Fs', ...
         'bits_per_sample', 'scale_factor', 'dict_custom', ...
         'nodes_left', 'nodes_right', 'nodes_sym', 'root');

    % 12) Print summary
    fprintf('\n===== Stage 1 Summary =====\n');
    fprintf('Entropy H (bits/symbol)            : %.4f\n', H);
    fprintf('Theoretical bits (H * N)           : %.0f\n', N_min_bits);
    fprintf('Original length (bits)             : %d\n', L_orig);
    fprintf('Huffman encoded length (bits)      : %d\n', L_huff);
    fprintf('Compression ratio (custom)         : %.4f\n', CR_huff);
    if ~isnan(L_builtin)
        fprintf('Compression ratio (built-in)       : %.4f\n', CR_builtin);
    end
    fprintf('Theoretical compression ratio K/H  : %.4f\n', CR_theory);
    fprintf('Reconstructed audio: stage1_reconstructed_custom.wav\n');
end
