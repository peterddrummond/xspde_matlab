function xwrite(errors,input,cdata,raw)
%   XWRITE(errors,input,data,raw) writes files produced by xsim.
%   Input:  'errors' vector,'input' parameter cell array,
%   'data' cell array, 'raw' trajectories cell array.
%   Output: HDF5 file or matlab file
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License.txt

[filename,hflag] = xwritename(input{1}.file);         %%get valid filename
if ~hflag                                             %%if file not blank
     save (filename,'errors','input','cdata','raw');  %%save data in Matlab file
     return
end                                                   %%end if file not blank
sequence = length(cdata);                             %%get sequence length
for s = 1:sequence                                    %%loop over sequence
    seq = sprintf('/data/sequence_%d',s);             %%name for sequence  
    graphs = size(cdata{s},4);                        %%get data length
    for g = 1:graphs                                  %%loop over graphs
        graphname = sprintf('/graph_%d',g);           %%define graph name
        dsname = [seq graphname];          
        h5create(filename, dsname, size(cdata{s}(:,:,:,g)));
        h5write(filename, dsname, cdata{s}(:,:,:,g));
    end                                               %%end loop over graphs
    h5writeatt(filename, seq, 'Graphs', graphs);      %%write graph number
end                                                   %%end loop over sequence
h5writeatt(filename, '/', 'Date', date());            %write attributes
h5writeatt(filename, '/', 'xSPDE_version', input{1}.version);
h5writeatt(filename, '/', 'Sequence', sequence);
for s = 1:sequence                                    %%loop over sequence
    xh5writecells(filename, [seq '/input'], input{s});%%input for sequence
end
end

function xh5writecells(filename, path, in)
%   XH5WRITECELLS(filename, inputname, in) writes HDF5 cell data.
%   Input:  file 'filename', data 'path', data source 'in'
%   Output: HDF5 file attributes including cell data
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License.txt
 
xh5writegroup(filename, path);
fields = fieldnames(in);
for i = 1:numel(fields)
    if iscell(in.(fields{i}))
        subpath = strcat(path,'/', fields{i});
        xh5writegroup(filename, subpath);
        acell = in.(fields{i});
        for j = 1:max(size(acell)) 
            xh5writeatt(filename, subpath, [fields{i} '_' int2str(j)], acell{j});
        end
    else
        xh5writeatt(filename, path, fields{i}, in.(fields{i}));
    end
end
end

function xh5writegroup(filename, path)
%   XH5WRITEGROUP(filename, inputname, in) creates empty HDF5 group.
%   Input:  file 'filename', data 'path'
%   Output: HDF5 file with new attribute group.
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License.txt
 
    plist = 'H5P_DEFAULT';
    fid = H5F.open(filename, 'H5F_ACC_RDWR', plist);
    gid = H5G.create(fid,path,plist,plist,plist);
    H5G.close(gid);
    H5F.close(fid);
end

function xh5writeatt(filename,location,attname,attvalue)
%   XH5WRITEATT(filename,location,attname,attvalue) 
%   Writes HDF5 data attributes produced by xsim.
%   Allows for functions by string conversion
%   Output: HDF5 file attribute 
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License.txt
 
if isa(attvalue,'function_handle')
    attvalue= ['function_' func2str(attvalue)];
end
h5writeatt(filename,location,attname,attvalue);
end

function [filename,hflag] = xwritename(filename)
%   [filename,hflag] = XWRITENAME(in_fname) 
%   Generates a valid filename for writing.
%   Returns hflag = 1 if file is an HDF5 file; 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt
 
hflag = 0;
[pathstr,name,ext] = fileparts(filename);
if strcmp(ext, '.h5')||strcmp(ext, '.mat')
    counter = 1;
    newname = name;
    while exist(fullfile(pathstr,[newname ext]), 'file')
        newname = sprintf('%s_%d', name, counter);
        counter = counter + 1;
    end
    filename = fullfile(pathstr,[newname ext]);
end
if strcmp(ext, '.h5')
    hflag = 1;
end
end