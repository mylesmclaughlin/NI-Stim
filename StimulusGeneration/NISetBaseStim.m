function NISetBaseStim(command)
global S

switch command
    case 'off'
        % delay normal stimulus
        S.stim.burstdelay = 0;
        % settings for base stimulus
        S.basestim.stim = 0;
    case 'frequency'
        % delay normal stimulus
        S.stim.burstdelay = 125;
        % settings for base stimulus
        S.basestim.stim = 1;
        S.basestim.amplitude = 1;
        S.basestim.frequency = 2;
        S.basestim.phase = 0;
        S.basestim.ampmoddepth = 0;
        S.basestim.ampmodfreq = 5;
        S.basestim.ampmodphase = 0.5;
        S.basestim.waveformindex = 1;
        S.basestim.phase1pulsewidth = 50;
        S.basestim.phase2pulsewidth = 50;
        S.basestim.phase1amp = 100;
        S.basestim.phase2amp = -100;
        S.basestim.phasegap = 0;
        S.basestim.waveformindex = 1;
    case 'ampmod'
        % settings for base stimulus
        S.basestim.stim = 1;
        S.basestim.amplitude = 1;
        S.basestim.frequency = 300;
        S.basestim.phase = 0;
        S.basestim.ampmoddepth = 100;
        S.basestim.ampmodfreq = 4;
        S.basestim.ampmodphase = 0.5;
        S.basestim.waveformindex = 1;
        S.basestim.phase1pulsewidth = 50;
        S.basestim.phase2pulsewidth = 50;
        S.basestim.phase1amp = 100;
        S.basestim.phase2amp = -100;
        S.basestim.phasegap = 0;
        S.basestim.waveformindex = 1;
        
        % delay normal stimulus
        S.stim.burstdelay = 125 - 1000*(1/300)/4;% delay normal stimulus
  
        
end



