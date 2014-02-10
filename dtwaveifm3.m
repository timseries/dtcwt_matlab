function X=dtwaveifm3(wlx,whx,nlevels)
    %does a 3-D DTCWT
     global xa ya h0o h1o g0o g1o h0a h1a g0a g1a
    
     szxa=size(ya{1}{1});%original size of xa
     if any(rem(szxa,2^nlevels))
        ext_mode = 8;    % 8 means the LoLoLo band at each level is to be appended to be a multiple of 8
     else
        ext_mode = 4;    % 4 means the LoLoLo band at each level is to be appended to be a multiple of 4
     end   
     xa=wlx;
     ya=whx;
         
     for k=[-nlevels:-1]; oct_cplx_xa(k); dtcwt3dC_xa(k,ext_mode); end   % Do nlevels levels of inverse transform
%      dtcwt3dCL_xa(-1,ext_mode);
    if nlevels==1
        rs=[2 1 1];
    else
        if nlevels==2
            rs=[2 2 1];
        else
            rs=[2 2 2];
        end
    end
     sxa=size(xa)/2; xa=xa(1:sxa(1),1:sxa(2),1:sxa(3));
     X=xa;
end