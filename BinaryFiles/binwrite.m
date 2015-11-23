function Done = binwrite(title,header,data,DataType,closeFile)

% BINWRITE - WRITE A BINARY FILE
%   DONE = BINWRITE(TITLE,HEADER,DATA,DATATYPE) - write DATA to a binary file called TITLE
%   with HEADER. DONE is 1 is writing was successful and 0 otherwise.
%
%   FID = BINWRITE(TITLE,HEADER,DATA,DATATYPE,CLOSEFILE) - giving the last
%   argument CLOSEFILE the value 0 means that the binary file is not closed
%   and the file ID is returned in FID. This option can facilitate
%   faster appending to the file. Appending to the file will also work when
%   the file is closed but may be slower.
%
%   For example,
%       binwrite('jump.bin','go get them',[1 2 3; 4 5 6; 7 8 9;],'double')
%
% See also binappend, binread and binfileinfo

% MMCL 28/03/08

% set params;
if nargin<4
    DataType = 'single';  %'double'; %'int16';
end
if nargin<5
    closeFile = 1;  %close file or leave open for faster appending
end

%append data params to header
lheader = length(header);
if lheader>999
    warning('Header must be less than 999 characters')
    Done = 0;
    return
end
dsize = size(data);
ncols = dsize(2);
if ncols>999
    warning('Data must have less than 999 coloums')
    Done = 0;
    return
end

% conver number of coloums, data type, subFileNum and length of header
MinInd = strfind(title,'-');
if isempty(MinInd)
    subFileNum = 0;
else
    DotInd = strfind(title,'.');
    subFileNum = str2num(title(MinInd+1:DotInd));
end
ncols = localconvert(ncols);
conDataType = localconvert(DataType);
postheader = ['/ncols=' ncols]; 
postpostheader = ['/dtype=' conDataType]; 
consubFileNum = localconvert(subFileNum);
postpostpostheader = ['/subfilenum=' consubFileNum];
lheader = localconvert(lheader);
preheader = ['nhead=' lheader '/'];

% make complete header
completeheader = [preheader header postheader postpostheader postpostpostheader '/'];

% update header length and remake header
lheader = length(completeheader);
lheader = localconvert(lheader);
preheader = ['nhead=' lheader '/'];
header = [preheader header postheader postpostheader postpostpostheader '/'];

%write file 
fid = fopen(title,'wb');  %open file in binary write mode
fwrite(fid, header, 'char');   %insert a header
fwrite(fid, data', DataType);   %write data to file
if closeFile
    fclose(fid);   %close file
    Done = 1;
else
    fclose(fid);   %close file and open in append mode
    fid = fopen(title,'a');
    Done = fid;
end

%------------------ Locals -----------------
function lstr = localconvert(lnum)

if ~ischar(lnum)
    lstr = num2str(lnum);
    while length(lstr)<3
        lstr = ['0' lstr];
    end  
    
%     if length(lstr)<3
%         if length(lstr)<2
%             lstr = ['00' lstr];
%         else
%             lstr = ['0' lstr];
%         end
%     end
else
    lstr = lnum;
    while length(lstr)<6
        lstr = [' ' lstr];
    end  
    
end

