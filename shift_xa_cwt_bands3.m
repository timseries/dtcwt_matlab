function shift_xa_cwt_bands3( hi_or_lo, level,motion_x,motion_y,motion_z,interpmethod )
% SHIFT_CWT_BANDS3 

% Modify the 28 subbands of CWT coefs in xa so as to be equivalent to
% shifting the source image by the vector field.
% The motion vector is described by motion_x, motion_y and motion_z.
% The motion vector is specified in units of voxel peroids at the current
% wavelet level.
% Note that xa contains all the 28 subbands and the LLL bands which is not
% converted into complex numbers. Complex bandpass interpoltion is used to
% process the 28 subbands.
% Huizhong Chen, Cambridge University, Oct 2008

global xa ya % xa is the LLL band, and ya is the highpass bands at each level

if nargin < 6, interpmethod = 'cubic'; end %default interpolation method is cubic
interpmethod_temp=interpmethod;
size_subband = size(ya{level}{1}); %sxa_subband is the size of 1 of the 8 subbands
s2 = size_subband(1)/2;
s1 = size_subband(2)/2;
s3 = size_subband(3)/2;

if level==1 && strcmp(hi_or_lo, 'hi')
    error('level 1 DTCWT contains no bandpass subbands');
end

% check if motion_x, motion_y and motion_z have the same dimension
if  isequal(size(motion_x), size(motion_y), size(motion_z)) ~= 1
    error('The dimensions of motion_x, motion_y and motion_z must agree');
end

if strcmp(hi_or_lo, 'hi')
    
    % check if the motion vectors have the same size as xa
    if isequal([s2 s1 s3], size(motion_x)) ~= 1
        error('subband must be of the same size as the motion vector');
    end

    % Set up matrix extension sizes and padding arrays at the specific level
    extx = ceil(max(abs(motion_x(:))));
    exty = ceil(max(abs(motion_y(:))));
    extz = ceil(max(abs(motion_z(:))));
    % define the padding matrices for x and y directions
    z1 = zeros(exty,s1+2*extx);
    z2 = zeros(s2,extx);
    % The padding in Z direction is done later

    % Create linear ramp matrices for phase wrapping
    thx_2d = ones(s2,1) * [1:s1] + extx; 
    thx = zeros(size(thx_2d,1),size(thx_2d,2),s3);
    for k=1:s3, thx(:,:,k)=thx_2d; end
    thy_2d = [1:s2]' * ones(1,s1) + exty; 
    thy = zeros(size(thy_2d,1),size(thy_2d,2),s3);
    for k=1:s3, thy(:,:,k)=thy_2d; end
    thz_2d = [1:s3]' * ones(1,s1) + extz; 
    thz = zeros(s2,s1,s3);
    for k=1:s3,
        for i=1:s2
        thz(i,:,k)=thz_2d(k,:);
        end
    end

    % Create matrices of interpolated point locations in the extended input
    % matrix.
    xs = motion_x + thx;
    ys = motion_y + thy;
    zs = motion_z + thz;

    % Set up expected rotation rates for each subband with freq offset.
    % To centre the passband of h on the scaling func and wavelet passbands.
    jw0 = -sqrt(-1) * pi/2.1;
    w = [jw0 jw0*3];


    w28 =  [w(2) w(1) w(1); -w(2) w(1) w(1); w(2) -w(1) w(1); -w(2) -w(1) w(1);  %HLL sub-band
            w(1) w(2) w(1); -w(1) w(2) w(1); w(1) -w(2) w(1); -w(1) -w(2) w(1);  %LHL sub-band
            w(2) w(2) w(1); -w(2) w(2) w(1); w(2) -w(2) w(1); -w(2) -w(2) w(1);  %HHL sub-band
            w(1) w(1) w(2); -w(1) w(1) w(2); w(1) -w(1) w(2); -w(1) -w(1) w(2);  %LLH sub-band
            w(2) w(1) w(2); -w(2) w(1) w(2); w(2) -w(1) w(2); -w(2) -w(1) w(2);  %HLH sub-band
            w(1) w(2) w(2); -w(1) w(2) w(2); w(1) -w(2) w(2); -w(1) -w(2) w(2);  %LHH sub-band
            w(2) w(2) w(2); -w(2) w(2) w(2); w(2) -w(2) w(2); -w(2) -w(2) w(2)]; %HHH sub-band

    % Get the subband which is denoted by index number 1 to 28
    % 1-4: HLL band; 5-8: LHL band; 9-12: HHL band; 13-16: LLH band; 
    % 17-20: HLH band; 21-24: LHH band; 25-28: HHH band
    for d = 1:28    
        x_extracted = extract_subband_xa( level, d ); % Extract the appropiate subband from ya
        x_extracted_converted_to_cpx = conv_cpx( x_extracted );
        x_extracted_e = zeros((2*size(z1,1)+size(x_extracted_converted_to_cpx,1)),(2*size(z2,2)+size(x_extracted_converted_to_cpx,2)),s3);
        % Having obtained the extracted_subband, unwrap the phases, and
        % extend borders in x and y directions.
        for k=1:s3
            %x_extracted_e(:,:,k) = [z1; z2 x_extracted_converted_to_cpx(:,:,k).*exp(-thx_2d*w28(d,1)-thy_2d*w28(d,2)-thz(:,:,k)*w28(d,3)) z2; z1];
            x_extracted_e(:,:,k) = [z1; z2 x_extracted_converted_to_cpx(:,:,k).*exp(-thx_2d*w28(d,1)-thy_2d*w28(d,2)-thz(:,:,k)*w28(d,3)) z2; z1];
        end

        % Now extend x_extracted_e in z direction
        size_x_extracted_e = size(x_extracted_e);
        z3 = zeros(size_x_extracted_e(1),size_x_extracted_e(2),extz);
        x_extracted_e = cat(3,z3,x_extracted_e,z3);  

        %must use linear interpolation if volume contains only 2 slices,
        %for example
        if min(size(x_extracted_e)) < 3
            interpmethod='linear';
        end
        % Interpolate x_extracted_e to the new points, specified in (xs,ys,zs).    
            x_extracted_i = interp3(x_extracted_e,xs,ys,zs,interpmethod);
        interpmethod=interpmethod_temp;  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rewrap the phases.
        for k=1:s3
             x_extracted_converted_to_cpx(:,:,k) = x_extracted_i(:,:,k).* exp(xs(:,:,k)*w28(d,1)+ys(:,:,k)*w28(d,2)+zs(:,:,k)*w28(d,3));
         end

        % Convert the complex numbers back to real and imag coefs, in order to
        % be stored back into xa
        x_extracted_converted_to_coefs = conv_coefs(x_extracted_converted_to_cpx);

        % Now update the highpass bands of ya by putting x_extracted at the correct
        % position in ya.
        update_subband_xa( x_extracted_converted_to_coefs, level, d );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LLL band interpolation - no need for any phase wrapping here.
    elseif strcmp(hi_or_lo, 'lo')
        s2 = size_subband(1);
        s1 = size_subband(2);
        s3 = size_subband(3);

        % Set up matrix extension sizes and padding arrays
        extx = ceil(max(abs(motion_x(:))));
        exty = ceil(max(abs(motion_y(:))));
        extz = ceil(max(abs(motion_z(:))));
        
        % Create linear ramp matrices for phase wrapping
        thx_2d = ones(s2,1) * [1:s1] + extx; for k=1:s3, thx(:,:,k)=thx_2d; end
        thy_2d = [1:s2]' * ones(1,s1) + exty; for k=1:s3, thy(:,:,k)=thy_2d; end
        thz_2d = [1:s3]' * ones(1,s1) + extz; 
        thz = zeros(s2,s1,s3);
        for k=1:s3,
           for i=1:s2
           thz(i,:,k)=thz_2d(k,:);
           end
        end

        % Create matrices of interpolated point locations in the extended input
        % matrix.
        xs = motion_x + thx;
        ys = motion_y + thy;
        zs = motion_z + thz;

        % Get the subband and extend its borders using symmetric extension.
        x_extracted_e = xa([[exty:-1:1] [1:s2] s2+1-[1:exty]],[[extx:-1:1] [1:s1] s1+1-[1:extx]],[[extz:-1:1] [1:s3] s3+1-[1:extz]]);
        
        if min(size(x_extracted_e)) < 3
            interpmethod='linear';
        end
        % Interpolate x_extracted_e to the new points, specified in (xs,ys,zs).
        x_extracted_i = interp3(x_extracted_e,xs,ys,zs,interpmethod);
        interpmethod=interpmethod_temp;
        % Now update the x matrix by putting x_extracted at the correct
        % position in xa.
        xa = x_extracted_i;
    else error('Input parameter hi_or_lo must be either hi or lo');
end

return