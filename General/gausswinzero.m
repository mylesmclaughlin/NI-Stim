function out = gausswinzero(samp)

out = gausswin(samp);
out = out-min(out);
out = out/max(out);

