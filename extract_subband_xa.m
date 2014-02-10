function  xa_extracted = extract_subband_xa( level, index_no )
%   xa_extracted = extract_subband_xa( level, index_no )
%   Extract the 1 of the 8 octals from xa in the specified level, 
%   e.g. subband = 'LLL', 'LHL', ect.
%   index_no specifies the which of the 4 subbands within the current octal
%   is selected. 
%   index_no = 1 => j1=+j, j2=+j, j3=+j;
%   index_no = 2 => j1=-j, j2=+j, j3=+j;
%   index_no = 3 => j1=+j, j2=-j, j3=+j;
%   index_no = 4 => j1=-j, j2=-j, j3=+j;
%   N.B. This version is only compatible with dtcwtDCL

global ya

if level<=0, error('level of DTCWT must be larger than 0'); end
if (index_no>28 || index_no<1)==1, error('index_number must be in the range of 1 to 28'); end

size_subband = size(ya{level}{1}) ; %size_subband is the size of 1 of the octal subbands

% Here we should add in a few lines to check if the subband has even size
if any(rem(size_subband,2)),
    error('Size of subband must be a multiple of 2');
end
s2 = size_subband(1)/2; s1 = size_subband(2)/2;

switch index_no
    case 1
        xa_extracted = ya{level}{1}(1:s2,1:s1,:);
    case 2
        xa_extracted = ya{level}{1}(1:s2,(s1+1):end,:);
    case 3
        xa_extracted = ya{level}{1}((s2+1):end,1:s1,:);
    case 4
        xa_extracted = ya{level}{1}((s2+1):end,(s1+1):end,:);
    case 5
        xa_extracted = ya{level}{2}(1:s2,1:s1,:);
    case 6
        xa_extracted = ya{level}{2}(1:s2,(s1+1):end,:);
    case 7
        xa_extracted = ya{level}{2}((s2+1):end,1:s1,:);
    case 8
        xa_extracted = ya{level}{2}((s2+1):end,(s1+1):end,:);
    case 9
        xa_extracted = ya{level}{3}(1:s2,1:s1,:);
    case 10
        xa_extracted = ya{level}{3}(1:s2,(s1+1):end,:);
    case 11
        xa_extracted = ya{level}{3}((s2+1):end,1:s1,:);
    case 12
        xa_extracted = ya{level}{3}((s2+1):end,(s1+1):end,:);
    case 13
        xa_extracted = ya{level}{4}(1:s2,1:s1,:);
    case 14
        xa_extracted = ya{level}{4}(1:s2,(s1+1):end,:);
    case 15
        xa_extracted = ya{level}{4}((s2+1):end,1:s1,:);
    case 16
        xa_extracted = ya{level}{4}((s2+1):end,(s1+1):end,:);
    case 17
        xa_extracted = ya{level}{5}(1:s2,1:s1,:);
    case 18
        xa_extracted = ya{level}{5}(1:s2,(s1+1):end,:);
    case 19
        xa_extracted = ya{level}{5}((s2+1):end,1:s1,:);
    case 20
        xa_extracted = ya{level}{5}((s2+1):end,(s1+1):end,:);
    case 21
        xa_extracted = ya{level}{6}(1:s2,1:s1,:);
    case 22
        xa_extracted = ya{level}{6}(1:s2,(s1+1):end,:);
    case 23
        xa_extracted = ya{level}{6}((s2+1):end,1:s1,:);
    case 24
        xa_extracted = ya{level}{6}((s2+1):end,(s1+1):end,:);
    case 25
        xa_extracted = ya{level}{7}(1:s2,1:s1,:);
    case 26
        xa_extracted = ya{level}{7}(1:s2,(s1+1):end,:);
    case 27
        xa_extracted = ya{level}{7}((s2+1):end,1:s1,:);
    case 28
        xa_extracted = ya{level}{7}((s2+1):end,(s1+1):end,:);
end

return;
