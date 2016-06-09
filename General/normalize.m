function normsig = normalize(sig)
normsig = sig-min(sig);
normsig = normsig/max(normsig);
normsig = normsig*2-1;

