function Yh_out=ri2c(Yh)
% function o2c(level)
% Convert the spatially-separated ouputs of oct_cplx_*() to a 4-D complex representation, with the 4th dimension indexing the subband.
% Maths annotations correspond to the definitions in Chen/Kinsgsbury "Efficient Registration of Nonrigid 3D Bodies"
    levs=length(Yh);
    subband_per_level=28; %4 directional subbands * (8-1) subbands at each level
    Yh_out=cell(levs,1);
%    Yh_z=cell(levs*subband_per_level,1);
    i=sqrt(-1);
    for level=1:levs
        Yhsz=size(Yh{level}{1})./[2 2 2];
        Yh_out{level}=zeros([Yhsz subband_per_level]);
        for j=1:subband_per_level
            index=(level-1)*subband_per_level+j;            
            index_band=floor((j-.01)/4)+1; %goes from 1 to 7
            Yh_temp=Yh{level}{index_band};
            wj{index}=zeros(Yhsz);
            endindex=size(Yh_temp)/2;
            switch rem(j,4)
              case 0 %$\Psi_4(x,y)$
                Yh_temp=Yh_temp((endindex(1)+1):1:end,(endindex(2)+1):1:end,1:2:end) + ...
                        i*Yh_temp((endindex(1)+1):1:end,endindex(2)+1:1:end,2:2:end);
              case 1 %$\Psi_1(x,y)$
                Yh_temp=Yh_temp(1:1:endindex(1),1:1:endindex(2),1:2:end) + ...
                        i*Yh_temp(1:1:endindex(1),1:1:endindex(2),2:2:end);
              case 2 %$\Psi_2(x,y)$
                Yh_temp=Yh_temp(1:1:endindex(1),(endindex(2)+1):1:end,1:2:end) + ...
                        i*Yh_temp(1:1:endindex(1),(endindex(2)+1):1:end,2:2:end);
              case 3 %$\Psi_3(x,y)$
                Yh_temp=Yh_temp((endindex(1)+1):1:end,1:1:endindex(2),1:2:end) + ...
                        i*Yh_temp((endindex(1)+1):1:end,1:1:endindex(2),2:2:end);
            end
            Yh_out{level}(:,:,:,j)=Yh_temp;
        end
    end    
end