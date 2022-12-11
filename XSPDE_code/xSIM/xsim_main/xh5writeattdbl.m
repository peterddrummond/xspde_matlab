function xh5writeattdbl(filename,path,attname,attvalue)
%   XH5WRITEATTDBL(filename,path,attname,attvalue)
%   Saves the attribute with name given by attname, value given by
%   attvalue to path given by path in HDF5 file given by filename
%   Attvalue is required to be of type double
%   The function call is delegated to h5writeatt (Matlab function).
%   If attvalue is a row-vector, we also add the suffix __rowvec to the name
%   to indicate this. The xread routine then takes this
%   into account and makes sure the read attribute becomes a row vector
%   again.
%   This is because h5writeatt will save row- or column-vectors as 1D data, thus
%   losing the information about if it was a row- or column-vector.
%   Called by: xwrite
%   Licensed by Peter D. Drummond, Simon Kiesewetter (2022) - see License.txt
    assert(isa(attvalue,'double'), 'xh5writeattdbl: attvalue is not of type double');
    if ismatrix(attvalue) && size(attvalue,1)==1 && size(attvalue,2)>1
        attname = [attname '__rowvec'];
    end
    h5writeatt(filename,path,attname,attvalue);
end
