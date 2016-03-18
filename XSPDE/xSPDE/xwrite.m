function input = xwrite(error,data,input,raw)
%   XWRITE(errors,input,data,raw) writes files produced by xsim.
%   Input:  'error' vector,'input' parameter cell array,
%   'data' cell array, 'raw' trajectories cell array.
%   Output: HDF5 file or matlab file
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License.txt

[filename,hflag] = xwritename(input{1}.file);         %%get valid filename
input{1}.file = filename;                             %%store filename 
switch hflag                                          %%check filename flag
case 0                                                %%if Matlab filename
  save (filename,'error','input','data','raw');       %%save  Matlab file
  return
case 1                                                %%if HDF5 filename
  sequence = length(data);                            %%get sequence length
  for s = 1:sequence                                  %%loop over sequence
    seq = sprintf('/data/sequence_%d',s);             %%name for sequence  
    graphs = size(data{s},4);                         %%get data length
    for g = 1:graphs                                  %%loop over graphs
        graphname = sprintf('/graph_%d',g);           %%define graph name
        dsname = [seq graphname];
        h5create(filename, dsname, size(data{s}(:,:,:,g)));
        h5write(filename, dsname, data{s}(:,:,:,g));
    end                                               %%end loop over graphs
    h5writeatt(filename, seq, 'Graphs', graphs);      %%write graph number
  end                                                 %%end loop over sequence
  h5writeatt(filename, '/', 'Date', date());          %write attributes
  h5writeatt(filename, '/', 'xSPDE_version', input{1}.version);
  h5writeatt(filename, '/', 'Sequence', sequence);
  for s = 1:sequence                                  %%loop over sequence
    seq = sprintf('/data/sequence_%d',s);             %%name for sequence  
    xh5writeattribute(filename, seq, 'input', input{s});%%input for sequence
  end
end
end

function xh5writeattribute(filename, path, attname, attvalue)
    if iscell(attvalue)                              %%is a cell array?
        subpath = strcat(path, '/', attname);        %%subpath name
        xh5writegroup(filename, subpath);            %%create subpath
        h5writeatt(filename, subpath, ...
            'XSPDE_iscell', '1');                    %%set cellarray flag
        for i = 1:max(size(attvalue))                %%loop over cells
           subsubpath = strcat(subpath, '/', ...     %%one sub-subgroup for              
               int2str(i));                          %%every cell
           xh5writegroup(filename, subsubpath);      %%create sub-subpath
           xh5writeattribute(filename, ...           %%recursively write cell
               subsubpath, 'value', attvalue{i});       
        end
        return;
    end
    if isstruct(attvalue)                            %%is a struct?
        subpath = strcat(path, '/', attname);        %%subpath name
        xh5writegroup(filename, subpath);            %%create subpath
        fields = fieldnames(attvalue);                   
        for i = 1:numel(fields)                      %%loop over fieldnames
            xh5writeattribute(filename, ...          %%recursively write fields
                subpath, fields{i}, attvalue.(fields{i}));
        end
        return;     
    end
    if isa(attvalue,'function_handle')               %is a function handle?
        attvalue= ['function_' func2str(attvalue)];  %store as string
    end
    h5writeatt(filename,path,attname,attvalue);      %write attribute here
end


%Version 1.03   xwrite prints error messages on filename errors

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

function [filename,hflag] = xwritename(filename)
%   [filename,hflag] = XWRITENAME(in_fname) 
%   Generates a valid filename for writing.
%   Returns hflag = 1 if file is an HDF5 file; 
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License.txt
 
hflag = 0;
[pathstr,name,ext] = fileparts(filename);
if strcmp(ext, '.h5')||strcmp(ext, '.mat')
    counter = 1;
    newname = name;
    while exist(fullfile(pathstr,[newname ext]), 'file')
        fprintf('Warning in xwritename: file %s exists\n',[newname ext]);
        newname = sprintf('%s_%d', name, counter);
        counter = counter + 1;
    end
    filename = fullfile(pathstr,[newname ext]);
    fprintf('Writing output to file %s\n',filename); %%output the filename
else
    fprintf('Error in xwritename: file %s has invalid type\n', filename);
    fprintf('Filename must end with .mat or .h5\n'); %%output error message
    hflag = -1;
end
if strcmp(ext, '.h5')
    hflag = 1;
end
end

%Version 1.03   xwritename prints error messages, returns an error flag