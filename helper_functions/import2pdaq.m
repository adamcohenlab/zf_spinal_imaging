function [data, srate] = import2pdaq(path,file,type) 
% this is meant to import voltage readings from a binary file. data are unscaled unsigned 16 or 32bit from labview
                                       % function expects
                                       %header with information about
                                       %number of channels and scaling
                                       %constants.
                                     % data are explicitly stored in
                                     % labview as big endian. file is
                                     % the filename. 
                                     %type= 'a', analog or 'c' counter

m = dir([path file]); % gives metadata, among other things the number of bytes in the file.
nbytes = m.bytes;

fid = fopen([path file]);
if(type == 'a')
h1 = fread(fid,2,'double','b'); %h1= first part of header; returns 2 constants: first one is number of scaling constants second one is number of channels). 
sc= fread(fid,h1(1)*h1(2),'double','b'); %sc=scalingconstants. one vector for all channels

sc=reshape(sc, h1(1), h1(2)); %scaling constants; 1 column for each channel. first value in each column is the sampling rate; 2nd to nth value are the scaling constants according to offset+x1*data+x2*data^2+x3*data^3 etc...
%%
h1(2)=1;                % header contains all AI info, only one channel per file saved
sc=sc(:,1);            %limit scaling constants to 1 channel: these two lines are specific for mode of saved analog channels, not generally correct
%%
spc=(nbytes-ftell(fid))/2/h1(2); %samples per channel. 32 bits=4 bytes; for 16 bit data divide by 2 instead of 4.ftell returns posiion in file which is where we are after reading out the header.

udat=fread(fid,[h1(2), spc],'uint16','b'); %read data from file into matrix with amount of columns equal to datapoints and amount or rows equal to amount of channels
sdat=udat-2^15;
fclose(fid);

temp=ones(spc,1)*sc(2,:);
for i=1:h1(1)-2

temp=temp+(ones(spc,1)*sc(i+2,:)).*(sdat'.^i);
end
data=temp;
srate=sc(1,1);
elseif(type == 'c')
   h1 = fread(fid,1,'double','b');
   srate=h1; 
   spc=(nbytes-ftell(fid))/2; %samples per channel. 32 bits=4 bytes; for 16 bit data divide by 2 instead of 4.ftell returns posiion in file which is where we are after reading out the header.

udat=fread(fid,[1, spc],'int16','b'); %read data from file into matrix with amount of columns equal to datapoints and amount or rows equal to amount of channels

fclose(fid);
   data=udat';
else
    'specify either analog or counter data'
end
end