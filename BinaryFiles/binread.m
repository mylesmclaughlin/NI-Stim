function D = binread(title,part,maxVarSize,warningon)

% BINREAD - READ BINARY FILE
%   D = BINREAD(TITLE) - read a binary file called TITLE and return data
%   in struct D. D.header contains the header of the file. D.data contains 
%   the data in a matlab matrix. D.part contains the part number of the file 
%
%   D = BINREAD(TITLE,PARTNUM) - requests a specific part of the file.
%   PARTNUM must be a positive integer. The PARTNUM will also be returned
%   in D.part. If the last part of the file is requested D.part will be
%   returend as zero.
%
%   D = BINREAD(TITLE,[IND1,IND2]) - requests a specific part of the file 
%   starting at IND1 and ending at IND2 for all coloums in D.data.
%   
%   D = BINREAD(TITLE,PARTNUM,MAXVARSIZE) - MAXVARSIZE sets the size of the
%   file parts. Default value is 1.5e7. Adjust if PC has a lot/little memory. 
%
%   For example,
%       D = binread('jump.bin',2)    
%
%   FILE FORMAT INFORMATION
%   The size of a binary file is limited 999GB. However, to keep files to a 
%   managable size data is stored in subfiles of 1GB. The subfile structure 
%   is handled internarly by binread and the user only ever needs to issue 
%   the command: binread('mybinaryfile.bin')
%
% See also binwrite, binappend and binfileinfo

% EVE	13/09/08 : modifications: 
%					* & -> &&; if nargin==3 && isempty(part) 
%					* dataL removed, direct concatination to reduce peak-memory.
%					* output is casted to stored datatype
% MMCL	28/03/08

% Check input arguments
if nargin==3 && isempty(part) 
    part = 1; % 
end
if nargin<2
    part = 1;
end
if nargin<3
    maxVarSize = 3e7; % will hit memory problems if variables exceed this size
end
if nargin<4
    warningon=1;
end
% get path
[pathstr,name,ext] =  fileparts(title);

% get header length - fisrt 9 chars
fid = fopen(title,'rb');  %open file
binnchars = fread(fid, 9, 'char');   %read length of header
fclose(fid);   %close file
nchars = char(binnchars(7:9)');
nchars = str2num(nchars);
nBytesHeader = nchars; % one byte per char ??

%read file
fid = fopen(title,'rb');   %open file
binheader = fread(fid, nchars, 'char');   %read in the header
fclose(fid); %close file

header = char(binheader');
ncolspos = strfind(header,'ncols=');
ncols = char(header(ncolspos+6:ncolspos+8));
ncols = str2num(ncols);
dtypepos = strfind(header,'dtype=');
DataType = header(dtypepos+6:dtypepos+11);
DataType = strtrim(DataType);
dtypepos = strfind(header,'subfilenum=');
if isempty(dtypepos)
    subFileNum = 0;
else
    subFileNum = char(header(dtypepos+11:dtypepos+13));
    subFileNum = str2num(subFileNum);
end
slash = strfind(header,'/');
header = header(slash(1)+1:slash(2)-1);

% get number of bytes per number stored in file
if strcmpi(DataType,'single')
    BytesPerNum = 4;
elseif strcmpi(DataType,'double')
    BytesPerNum = 8;    
elseif strcmpi(DataType,'int32')
    BytesPerNum = 4;
else
    if warningon;  disp(['Unrecognized data type: ' DataType]); end
    D = 0;
    return
end

% get subfile info 
F = binfileinfo(title);
totalDataSize = 0;
for n = 1:length(F)
    F(n).nRows =   (F(n).bytes - nBytesHeader)/BytesPerNum/ncols;
    totalDataSize = totalDataSize + (F(n).bytes - nBytesHeader);
end
totalnRows = totalDataSize/BytesPerNum/ncols;

% determine the starting and ending row of the requested data
if length(part) == 1
    rowStart = (part-1) * floor(maxVarSize/ncols) + 1;
    rowEnd = rowStart + floor(maxVarSize/ncols) - 1;
elseif length(part) == 2
    rowStart = part(1);
    rowEnd = part(2);
else
     if warningon; disp(['Part input not correctly formatted: ' num2str(part)]); end
    D = 0;
    return
end

% check maxVarSize is not exceeded
if (rowEnd - rowStart)*ncols>maxVarSize
     if warningon; disp(['Data requested exceeds maxVarSize: ' num2str(maxVarSize)]); end
    D = 0;
    return
end

% Does the data request hit the end of the file?
if rowEnd>totalnRows 
   rowEnd = totalnRows;
   if length(part) == 1
        part = 0;
   elseif length(part) == 2
       part(2) = rowEnd;
   end
end

% which subfiles will are these rows contained in? 
lastrow = 0;
for n = 1:length(F)
    lastrow = lastrow + F(n).nRows;
    if rowStart<=lastrow
        FileStart = n;
        break
    end
end
if ~exist('FileStart','var')
     if warningon; disp('Data requested exceeds filesize'); end
    D = 0;
    return
end
lastrow = 0;
for n = 1:length(F)
    lastrow = lastrow + F(n).nRows;
    if rowEnd<=lastrow
        FileEnd = n;
        break
    end
end

% Start reading in the data
data = [];


rowcount = 0;
if FileStart > 1
    for n = 1:FileStart-1
        rowcount = rowcount + F(n).nRows;
    end
end
for n = FileStart:FileEnd
    if n == FileStart % determine start row for this subfile
        LrowStart = rowStart - rowcount;
    else
        LrowStart = 1;
    end
    if n == FileEnd % determine end row for this subfile
        LrowEnd = rowEnd - rowcount;
    else
        LrowEnd = F(n).nRows;
    end
    
     data = cat(1, data, (localbinread([pathstr '\' F(n).name],LrowStart,LrowEnd,ncols,nBytesHeader,DataType,BytesPerNum))); % get the data from the correct subfile
  
    % keep count of the number of rows already read
    rowcount = rowcount + F(n).nRows;
end

%make return struct
D.header		= header;
D.data			= data;
D.part			= part;
D.subFileNum	= [FileStart:FileEnd]-1;
D.maxVarSize	= maxVarSize;

% ----------------------- Locals ----------------------
function data = localbinread(title,rowStart,rowEnd,ncols,nBytesHeader,DataType,BytesPerNum);

% open file
fid = fopen(title,'rb');   

% define offset and number of rows
offset	= (rowStart-1)*ncols*BytesPerNum + nBytesHeader;
nrows	= (rowEnd-rowStart)+1;

% Get the data
fseek(fid, offset, 0);						% set file position indicator
data = fread(fid, [ncols nrows], DataType); %read in the data
data = data';								%must transpose data after reading in
data = cast(data,DataType);					%cast to real class
fclose(fid);								%close file
