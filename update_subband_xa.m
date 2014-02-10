function  update_subband_xa( update_matrix, level, index_no )
% UPDATE_SUBBAND
% This function is to update the apporipriate section of the subband in xa
% with update_matrix
global ya

if level<=0, error('level of DTCWT must be larger than 0'); end
if (index_no>28 || index_no<1)==1, error('index_number within the current subband must be in the range of 1 to 4'); end

size_subband = size(ya{level}{1}); %size_subband is the size of 1 of the octal subbands
s2 = size_subband(1)/2; s1 = size_subband(2)/2; %s3 = size_subband(3)/2;
% Here we should add in a few lines to check if the subband has even size
if any(rem(size_subband,2)),
    error('Size of subband must be a multiple of 2');
end

switch index_no
    case 1
        ya{level}{1}(1:s2,1:s1,:) = update_matrix;
    case 2
        ya{level}{1}(1:s2,(s1+1):end,:) = update_matrix;
    case 3
        ya{level}{1}((s2+1):end,1:s1,:) = update_matrix;
    case 4
        ya{level}{1}((s2+1):end,(s1+1):end,:) = update_matrix;
    case 5
        ya{level}{2}(1:s2,1:s1,:) = update_matrix;
    case 6
        ya{level}{2}(1:s2,(s1+1):end,:) = update_matrix;
    case 7
        ya{level}{2}((s2+1):end,1:s1,:) = update_matrix;
    case 8
        ya{level}{2}((s2+1):end,(s1+1):end,:) = update_matrix;
    case 9
        ya{level}{3}(1:s2,1:s1,:) = update_matrix;
    case 10
        ya{level}{3}(1:s2,(s1+1):end,:) = update_matrix;
    case 11
        ya{level}{3}((s2+1):end,1:s1,:) = update_matrix;
    case 12
        ya{level}{3}((s2+1):end,(s1+1):end,:) = update_matrix;
    case 13
        ya{level}{4}(1:s2,1:s1,:) = update_matrix;
    case 14
        ya{level}{4}(1:s2,(s1+1):end,:) = update_matrix;
    case 15
        ya{level}{4}((s2+1):end,1:s1,:) = update_matrix;
    case 16
        ya{level}{4}((s2+1):end,(s1+1):end,:) = update_matrix;
    case 17
        ya{level}{5}(1:s2,1:s1,:) = update_matrix;
    case 18
        ya{level}{5}(1:s2,(s1+1):end,:) = update_matrix;
    case 19
        ya{level}{5}((s2+1):end,1:s1,:) = update_matrix;
    case 20
        ya{level}{5}((s2+1):end,(s1+1):end,:) = update_matrix;
    case 21
        ya{level}{6}(1:s2,1:s1,:) = update_matrix;
    case 22
        ya{level}{6}(1:s2,(s1+1):end,:) = update_matrix;
    case 23
        ya{level}{6}((s2+1):end,1:s1,:) = update_matrix;
    case 24
        ya{level}{6}((s2+1):end,(s1+1):end,:) = update_matrix;
    case 25
        ya{level}{7}(1:s2,1:s1,:) = update_matrix;
    case 26
        ya{level}{7}(1:s2,(s1+1):end,:) = update_matrix;
    case 27
        ya{level}{7}((s2+1):end,1:s1,:) = update_matrix;
    case 28
        ya{level}{7}((s2+1):end,(s1+1):end,:) = update_matrix;
end

return;