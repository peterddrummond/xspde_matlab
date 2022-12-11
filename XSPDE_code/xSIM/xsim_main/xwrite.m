function xwrite(error,data,input,raw)
%   XWRITE(error,data,raw,input) writes files produced by xsim.
%   Input:  error vector, 'data' cell array,'raw' trajectories cell array,
%   'input' parameter cell array.
%   Output: HDF5 file or matlab file
%   Note: this version does not save the error data in HDF5!
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License

hflag = xwriteprep(input{1}.file,input{1});           %%make copy if needed
filename = input{1}.file;                             %%store filename 
switch hflag                                          %%check filename flag
case 0                                                %%if Matlab filename
  save (filename,'error','data','input','raw');       %%save  Matlab file
  return
case 1                                                %%if HDF5 filename
  sequence = length(data);                            %%get sequence length
  for s = 1:sequence                                  %%loop over sequence
    seq = sprintf('/data/sequence_%d',s);             %%name for sequence  
    graphs = size(data{s},2);                         %%get data length
    for g = 1:graphs                                  %%loop over graphs
        graphname = sprintf('/graph_%d',g);           %%define graph name
        dsname = [seq graphname];
        h5create(filename, dsname, size(data{s}{g}));
        h5write(filename, dsname, real(data{s}{g}));  %%store real part 
    end                                               %%end loop on graphs
    h5writeatt(filename, seq, 'Graphs', graphs);      %%write graph number
  end                                                 %%end sequence loop 
  h5writeatt(filename, '/', 'Date', date());          %write attributes
  h5writeatt(filename, '/', 'xSPDE_version',  input{1}.version);
  h5writeatt(filename, '/', 'xGRAPH_version', input{1}.version);
  h5writeatt(filename, '/', 'Sequence', sequence);
  for s = 1:sequence                                  %%loop over sequence
   seq = sprintf('/data/sequence_%d',s);              %%name for sequence  
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
           xh5writeattribute(filename, ...           %%recursive write cell
               subsubpath, 'value', attvalue{i});       
        end
        return;
    end
    if isstruct(attvalue)                            %%is a struct?
        subpath = strcat(path, '/', attname);        %%subpath name
        xh5writegroup(filename, subpath);            %%create subpath
        fields = fieldnames(attvalue);                   
        for i = 1:numel(fields)                      %%loop over fieldnames
            xh5writeattribute(filename, ...          %%recursive write fields
                subpath, fields{i}, attvalue.(fields{i}));
        end
        return;     
    end
    if isa(attvalue,'function_handle')               %%is a function handle?
        attvalue= ['function_' func2str(attvalue)];  %%store as string
    end
    if isa(attvalue,'double')
        if ~isreal(attvalue)                         %%if it is a complex array
           xh5writeattdbl(filename,path,attname,real(attvalue));%%write real part
           xh5writeattdbl(filename,path,strcat(attname,'__imag'),imag(attvalue));%%write imaginary part
        else
            xh5writeattdbl(filename,path,attname,attvalue);      %%write real-valued array here
        end
    else  %%if attvalue is neither double nor cell array or struct
        h5writeatt(filename,path,attname,attvalue);
    end
end

%Version 1.03   xwrite prints error messages

function xh5writegroup(filename, path)
%   XH5WRITEGROUP(filename, inputname, in) creates empty HDF5 group.
%   Input:  file 'filename', data 'path'
%   Output: HDF5 file with new attribute group.
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015): see License
 
    plist = 'H5P_DEFAULT';
    fapl  = H5P.create('H5P_FILE_ACCESS');
    cfapl = onCleanup(@()H5P.close(fapl));
               %%this is important for large attributes (>64KB)
    H5P.set_libver_bounds(fapl,'H5F_LIBVER_18','H5F_LIBVER_LATEST');                                                    
    fid   = H5F.open(filename, 'H5F_ACC_RDWR', fapl);
    cf    = onCleanup(@()H5F.close(fid));
    gid   = H5G.create(fid,path,plist,plist,plist);
    gf    = onCleanup(@()H5G.close(gid)); 
end

function [hflag] = xwriteprep(filename,p)
%   [hflag] = XWRITEPREP(in_fname) 
%   Checks if filename exists and makes a copy if needed.
%   Returns hflag = 1 if file is an HDF5 file; 
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License
 
hflag = 0;
[pathstr,name,ext] = fileparts(filename);
if strcmp(ext, '.h5')||strcmp(ext, '.mat')
    if exist(filename, 'file')
       newname = [name '_1'];
       newfname = fullfile(pathstr,[newname ext]);
       xpr(0,p,'Warning in xwritename: file %s exists, copying to %s\n',...
           [name ext], [newname ext]);
       movefile(filename, newfname);
    end
    xpr(1,p,'Writing output to file %s\n',filename); %%output the filename
else
    xpr(-1,p,'Error in xwritename: file %s has invalid type\n', filename);
    xpr(-1,p,'Filename must end with .mat or .h5\n'); %%output error message
    hflag = -1;
end
if strcmp(ext, '.h5')
    hflag = 1;
end
end
