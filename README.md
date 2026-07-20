# Audio Source and Channel Coding

ECNG-4302 Project 2 — MATLAB simulation of an audio transmission system over a binary symmetric channel (BSC).

## Pipeline

1. **Source coding** — Quantize audio to 6-bit symbols, apply custom Huffman coding (tree built from symbol PMF), verify lossless reconstruction.
2. **Channel coding** — Transmit Huffman-compressed bits over a BSC with two schemes:
   - *Uncoded*: Huffman only, tested at flip probabilities 0–0.02.
   - *Coded*: Rate-1/2 convolutional code (constraint length 3, generators [7 5]) with hard-decision Viterbi decoding, tested at flip probabilities 0–0.05.

## Usage

Open in MATLAB, then:

```matlab
>> main
```

Select mode 1 (source coding only) or mode 2 (full pipeline). Output WAV files are saved to the current directory.

## Files

| File | Description |
|---|---|
| `main.m` | Mode selection entry point |
| `source_coding.m` | Huffman source coding stage |
| `channel_coding.m` | BSC + convolutional coding stage |
| `sample.wav` | Input audio example |
