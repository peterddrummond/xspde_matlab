function [hflag,input,data,raw] =  xread(gr_input)
%   [hflag,input,cdata,raw]  = XREAD(input) reads xSPDE data files. 
%   Input:  'filename' from input  
%   If input is cell array, reads data in file 'input{1}.file'. 
%   If input is string: 'fname', reads input and data from file 'fname'. 
%   Output: 'input' parameter cell array, data cell array 'cdata'.
%   Optional: 'raw' trajectories from a file. 
%   MIT licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License.txt

raw={};                                        %%initialize raw data
data = {};                                     %%set data to cell
saved = 0;                                     %%saved  = 0 for new input 
if ischar(gr_input)                            %%if input is character
        fname=gr_input;                        %%store label
        saved = 1;                             %%saved  = 1 for saved input 
else                                           %%end check if data is label
        input = xmakecell(gr_input);
        fname = input{1}.file;                 %%store label
end                                            %%end if data is character
hflag = xreadname(fname);                      %%check filename is OK
if hflag == -1                                 %%if filename not OK 
    return;                                    %%return to xgraph 
end                                            %%end if filename 
if hflag == 0                                  %%if filename  Matlab type 
    load(fname);                               %%load file data 
    try                                        %%test Matlab filename
        v = input{1}.version;                  %%get input version name
        fprintf('Reading Matlab file %s generated by %s\n',fname,v);
    catch                                      %%error in file filename
        fprintf('Warning in xread: file %s is not an xSPDE file\n', fname);
    end
    if saved == 0                              %% if new input required
        input = gr_input;                      %%use graphics input 
    end                                        %%end if  saved flag = 0
    return;                                    %%return to xgraph
end                                            %%end test Matlab filename
try                                            %%test HDF5 filename
    v = h5readatt(fname,'/','xSPDE_version');  %%get version name
    fprintf('Reading HDF5 file %s generated by %s\n', fname,v);
catch                                          %%error in file filename
    fprintf('Warning in xread: file %s is not an xSPDE file\n', fname);
end                                            %%end test HDF5 filename
sequence = h5readatt(fname,'/','Sequence');    %%define sequence number
for s = 1:sequence                             %%loop over sequence
    in = struct;                               %%define structure
    seq = sprintf('/data/sequence_%d',s);      %%name for sequence  
    graphs = h5readatt(fname, seq, 'Graphs');  %%define graph number
    for g = 1:graphs                           %%loop over graphs
        graphname = sprintf('/graph_%d',g);
        dsname = [seq graphname]; 
        data{s}(:,:,:,g)= h5read(fname,dsname);
    end                                        %%end loop graphs
    in=xh5attributes(fname, [seq '/input']);    
    input{s} = in;
end                                            %%end loop sequence
    if saved == 0                              %% if saved flag = 0
        input = gr_input;                      %%use graphics input 
    end                                        %%end if  saved flag = 0
end

function [res] = xh5attributes(filename, path)
    info = h5info(filename, path);
    num_groups = max(size(info.Groups));            %%number of subgroups
    num_attributes = max(size(info.Attributes));    %%number of attributes
    is_cell = 0;                                    %%cell array flag
    for i=1:num_attributes
       if strcmp(info.Attributes(i).Name, ...       %%if cell flag is set
               'XSPDE_iscell') 
           is_cell = 1;                             %%mark as cell array
       end
    end
    for i=1:num_groups
       if all(isstrprop(info.Groups(i).Name, ...    %%if have purely numeric
               'digit'))                            %%named subgroup
          is_cell = 1;                              %%mark as cell array
       end
    end
    if (~is_cell)                                   %%if not a cell array
        for i = 1:num_attributes                    %%loop over attributes
            attname = info.Attributes(i).Name;
            attvalue = info.Attributes(i).Value;
            if ~ischar(attvalue)
                attvalue = attvalue';
            else
                len = length (attvalue);            %%if is function handle
                if (len > 9) && ...                 %%read function as string
                        strcmp('function_', attvalue(1:9))
                    attvalue = str2func(attvalue(10:len));
                end
            end
        res.(attname) = attvalue;
        end
        for i = 1:num_groups                         %%loop over subgroups
            group_name = info.Groups(i).Name;
            tmp = xh5attributes(filename, group_name);%recursively read subgroup           
            res.(xcut_off_path(group_name)) = tmp;
        end
    else                                            %%if is a cell array
        res = {};
        for i=1:num_groups                          %%loop over subgroups
            subpath = strcat(path, '/', int2str(i));
            tmp = xh5attributes(filename, subpath); %%recursively read subgroup
            res{i} = tmp.value;
        end
    end            
end

%Version 1.03   xread returns an error flag on read error

function [res] = xcut_off_path(str)
%   XCUT_OFF_PATH(str) removes path data from an HDF5 group name
%   Output: HDF5 group name without path data. 
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License.txt
    C = strsplit(str,'/');
    res = C{size(C,2)};
end
    
function hflag = xreadname(filename)
%   [filename,hflag] = XREADNAME(in_fname) 
%   Tests for a valid input filename for reading.
%   Returns hflag = 1 if file is an HDF5 file, '.h5';
%   Returns hflag = 0 if file is a matlab file, '.mat';
%   Returns hflag = -1 if file does not exist or invalid.
%   Licensed by Peter D. Drummond & Simon Kiesewetter (2015) - see License.txt
 
hflag = -1;
[~,~,ext] = fileparts(filename);
if strcmp(ext, '.h5')
    hflag = 1;
end 
if strcmp(ext, '.mat')
    hflag = 0;
end 
if hflag == -1
    fprintf('Error in xreadname: file %s has invalid type\n', filename);
    fprintf('Filename must end with .mat or .h5\n'); %%output error message
else
  if ~exist(filename, 'file')
    fprintf('Error in xreadname: file %s doesn''t exist\n', filename);
    hflag = -1;
  end
end
end

%Version 1.03   xreadname has improved error message printing

