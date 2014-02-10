function [wlx whx]=dtwavexfm3(X,nlevels)
        %does a 3-D DTCWT
         global xa ya h0o h1o g0o g1o h0a h1a g0a g1a
        %  load3D_xa('sphere'); 
         xa=X;
         ya={};
         if any(rem(size(xa),2^nlevels))
            ext_mode = 8;    % 8 means the LoLoLo band at each level is to be appended to be a multiple of 8
         else
            ext_mode = 4;    % 4 means the LoLoLo band at each level is to be appended to be a multiple of 4
         end

        for k=1:3; xa=cat(k,xa,0*xa); end % Load a dataset and double its dimensions
%         dtcwt3dCL_xa(1,ext_mode);
        for k=1:nlevels; dtcwt3dC_xa(k,ext_mode); oct_cplx_xa(k); end   % Do nlevels levels of forward transform
        wlx=xa;
        whx=ya;
end