function D = binfileinfo(title)

% BINFILEINFO List info on subfile structure of bin file
%   F = BINFILEINFO(TITLE) returns a struct F containing info on all
%   subfiles related to the bin file TITLE. Each element of D contains info
%   on one subfile.
%   
%   F(n).name  -- subfilename
%   F(n).date  -- subfile modification date
%   F(n).bytes -- number of bytes allocated to the subfile
%
%   Example - 
%       F = binfileinfo('jump.bin')
%
%   % See also binwrite, binappend and binread

%   MMCL 10/07/2008

% check for existing files and subfiles 
if exist(title,'file')
    [pathstr,name,ext] =  fileparts(title);
    TrySubFileNum = '';
    n = 1;
    while exist(fullfile(pathstr, [name TrySubFileNum ext]),'file')
        SubFileNum = TrySubFileNum;
        title = fullfile(pathstr, [name SubFileNum ext]); % arry of subfiles
        D(n) = dir(title);
        n = n+1;
        if isempty(SubFileNum)
            TrySubFileNum = '-1';
        else
            NumTrySubFileNum = str2num(SubFileNum(2:end)) + 1;
            TrySubFileNum = ['-' num2str(NumTrySubFileNum)];
        end
    end
else
    warning([title ' does not exist.'])
    Done = 0;
    return
end

D = rmfield(D,'isdir');
