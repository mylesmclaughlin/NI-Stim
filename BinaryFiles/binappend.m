function Done = binappend(title,data,maxFileSize,fileID,DataType)

% BINWRITE - APPEND BINARY FILE
%   D = BINAPPEND(TITLE,DATA) - append DATA to a binary file called TITLE.
%   D returns the  number of the subfile to which the data was appended.
%   If appending was unsuccessful D is -1.
%
%   FILE FORMAT INFORMATION
%   The size of a binary file is limited 999GB. However, to keep files to a
%   managable size data is stored in subfiles of 1GB. If appending data to a
%   file (eg. mybinaryfile.bin) would make this file larger than 1GB then a
%   new binary subfile will be opened with a suffix -1 (eg. mybinaryfile-1.bin).
%   The subfile structure is handled internarly by binappend and the user only
%   ever needs to issue the command: binappend('mybinaryfile.bin',data)
%
%   BINAPPEND(TITLE,DATA,MAXFILESIZE) - specifies the MAXFILESIZE of the
%   subfiles in bytes. Defaulf size is 1GB
%
%   FID = BINAPPEND(TITLE,DATA,[],FID,DATATYPE) - if BINWRITE did not close the file
%   the FID can be given to facilitate faster appending to the file. The
%   FID is returned and is updated for each new subfile. DATATYPE must also
%   be given ('single','double','int32')
%
%   For example,
%       binappend('jump.bin',[1 2 3; 4 5 6; 7 8 9]);
%
% See also binwrite binread and binfileinfo

% MMCL 28/03/08

if nargin<3
    maxFileSize = 1e9; % 1GB : easier to work with files smaller the 1GB
end
if isempty(maxFileSize)
    maxFileSize = 1e9; % 1GB : easier to work with files smaller the 1GB
end
if nargin<4
    fileID = []; % fid to allow faster appending
end

% check for existing files and files with subnames
if exist(title,'file')
    [pathstr,name,ext] =  fileparts(title);
    TrySubFileNum = '';
    while exist(fullfile(pathstr,[name TrySubFileNum ext]),'file')
        SubFileNum = TrySubFileNum;
        if isempty(SubFileNum)
            TrySubFileNum = '-1';
        else
            NumTrySubFileNum = str2num(SubFileNum(2:end)) + 1;
            TrySubFileNum = ['-' num2str(NumTrySubFileNum)];
        end
    end
else
    warning([title ' does not exist. Try using binwrite'])
    Done = -1;
    return
end
% New current title
title = fullfile(pathstr,[name SubFileNum ext]);

% get header length - fisrt 9 chars
if isempty(fileID)
    fid = fopen(title,'rb');   %open file
    binnchars = fread(fid, 9, 'char');   %read length of header
    fclose(fid);   %close file
    nchars = char(binnchars(7:9)');
    nchars = str2num(nchars);
    % get number of coloums in bin file and data type
    fid = fopen(title,'rb');  %open file
    binheader = fread(fid, nchars, 'char');   %read in the header
    fclose(fid);   %close file

    header = char(binheader');
    ncolspos = strfind(header,'ncols=');
    ncols = char(header(ncolspos+6:ncolspos+8));
    ncols = str2num(ncols);
    dtypepos = strfind(header,'dtype=');
    DataType = header(dtypepos+6:dtypepos+11);
    DataType = strtrim(DataType);
    slash = strfind(header,'/');
    header = header(slash(1)+1:slash(2)-1);

    % Check data size matches data in file
    dsize = size(data);
    ncolsindata = dsize(2);
    if ncols~=ncolsindata
        warning(['Number of coloums in DATA does not match number of colums in ' title])
        Done = -1;
        return
    end
else
    fid = fileID;
    ss = size(data);
    ncols = ss(2);
end

% check the type of data and size of file
if strcmpi(DataType,'single')
    BytesPerNum = 4;
elseif strcmpi(DataType,'double')
    BytesPerNum = 8;
elseif strcmpi(DataType,'int32')
    BytesPerNum = 4;
else
    disp(['Unrecognized data type: ' DataType]);
    Done = -1;
    return
end

FileInfo = dir(title);
FileSize = FileInfo.bytes;
DataSize = length(data(:))*BytesPerNum;
if FileSize+DataSize<maxFileSize % just append to current file
    data1 = data;
    data2 = [];
elseif FileSize+DataSize>maxFileSize % split data over this current file and a new file
    Bytes1 = maxFileSize - FileSize;
    Nums1 = Bytes1/BytesPerNum;
    Rows1 = round(Nums1/ncols);
    if Rows1==0
        data1 = [];
    else
        data1 = data(1:Rows1,:);
    end
    data2 = data(Rows1+1:end,:);
end

% append data
if ~isempty(data1)
    %append file to this file
    if isempty(fileID)
        fid = fopen(title,'a');  %open file in binary write mode
    end
    
    fwrite(fid, data1', DataType);   %write data to file    
    if isempty(fileID)
        fclose(fid);   %close file
        if isempty(SubFileNum)
            NumSubFileNum = 0;
        else
            NumSubFileNum = str2num(SubFileNum(2:end));
        end
        Done = NumSubFileNum;
    else
        Done = fid;
    end
end

if ~isempty(data2)
    
    if ~isempty(fileID)% get header
        D = binread(title,[1 2]);
        header = D.header;
    end
    
    % wrtie data to new file
    if isempty(SubFileNum)
        NewSubFileNum = '-1';
        NumSubFileNum = 1;
    else
        NumSubFileNum = str2num(SubFileNum(2:end)) + 1;
        NewSubFileNum = ['-' num2str(NumSubFileNum)];
    end
    title = fullfile(pathstr,[name NewSubFileNum ext]);

    if isempty(fileID)
        Done = binwrite(title,header,data2,DataType);
        if Done == 1
            Done = NumSubFileNum;
        end
    else
        Done = binwrite(title,header,data2,DataType,0);
        fclose(fid);
    end
end
