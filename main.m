%% Audio Source and Channel Coding
% ECNG-4302 Project 2 - Full pipeline (Huffman + BSC + Convolutional coding)

clc; clear; close all;

fprintf('Simulation of Source and Channel Coding for an Audio System over a BSC\n');
fprintf('====================================================================\n');
fprintf('1) Stage 1 only: Huffman source coding\n');
fprintf('2) Full pipeline: Stage 1 + Stage 2\n');
fprintf('====================================================================\n');

mode = [];
while isempty(mode) || ~ismember(mode,[1 2])
    mode = input('Select mode (1 or 2): ');
end

if mode == 1
    fprintf('\n>> Running Stage 1 (Huffman source coding)...\n');
    source_coding();
    fprintf('>> Stage 1 completed.\n');
else
    fprintf('\n>> Running Stage 1 (Huffman source coding)...\n');
    source_coding();
    fprintf('>> Running Stage 2 (BSC + convolutional coding)...\n');
    channel_coding();
    fprintf('>> Full pipeline completed.\n');
end

fprintf('====================================================================\n');
fprintf('Simulation finished.\n');
fprintf('====================================================================\n');
