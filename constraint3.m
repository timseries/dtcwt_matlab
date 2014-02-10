function cvecs = constraint3(level, w28, mode, it)
% function cvecs = constraint3(level, w28, mode, it)

% Generate the constraint vectors cvecs from the reference xb and previous
% xa subband coefs.

% One set of 28 constaints is produced for each 2*2*2 group of subband
% coefs.

% w28 is the expected phase shift per sample for each subband which allows
% proper interpretation of measured phase shifts mod 2*pi.

% The units of shift are assumed to be (half) the image size and not the
% voxel spacing (since this changes from level to level).

% cvecs{d} is a 4-D array in (y,x,z,cvec element), where d is the subband number (1 to 28).
% constraints3 uses adjacent 2*2*2 groups of coefs so that size(cvecs,1:3)
% is just one less than size of the subband(28).

% The value of mode(1) indicates the number of global motion iterations 
% mode(2) = 1 removes the constraints for the outer edge of cvecs
% mode(2) = 2 removes the constraints for the outer edge of cvecs on only the
% subbands which a re not normal to the edge of the graph.

global ya

% if isequal(size(xa),size(xb),size(ya),size(yb))~=1
%     error('The input images must have the same size');
% end

sxa_subband = size(ya{level}{1}); %sxa_subband is the size of 1 of the 8 subbands
s2 = sxa_subband(1)/2;
s1 = sxa_subband(2)/2;
s3 = sxa_subband(3);

% CHECK HERE, WE ACTUALLY NEED TO MAKE SURE THE SUBBAND SIZE IS A MULTIPLE
% OF 2.
% if any(rem([s2 s1 s3],2))
%     error('subbands must be of even size for CONSTRAINT3.M');
% end

% sy, sx and sz specify the locality size of cvecs, indicating rows,
% columns and frames respectively
sy = s2-1;
sx = s1-1;
sz = s3/2-1;

ty = 1:sy;
tx = 1:sx;
tz = 1:sz;


for d = 1:28
    xa_extracted = extract_subband_xa( level, d ); % extract the subband
    prev_complex = conv_cpx(xa_extracted); % convert the DTCWT coefs to complex numbered coefs. As a result the number of frams is reduced by a factor of 2.
    xb_extracted = extract_subband_xb( level, d );
    ref_complex = conv_cpx(xb_extracted);  % convert the DTCWT coefs to complex numbered coefs. As a result the number of frams is reduced by a factor of 2.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each subband, measure horizontal phase gradients by taking the angle of summed
    % conjugate products across horizontal pairs.
    prevcp = prev_complex(:,tx+1,:) .* conj(prev_complex(:,tx,:));
    refcp = ref_complex(:,tx+1,:) .* conj(ref_complex(:,tx,:));
    
    % The next line calculates the horiztonal phase gradients in subband d.
    hcp = refcp(ty,:,tz) + refcp(ty+1,:,tz) + refcp(ty,:,tz+1) + refcp(ty+1,:,tz+1)...
        + prevcp(ty,:,tz) + prevcp(ty+1,:,tz) + prevcp(ty,:,tz+1) + prevcp(ty+1,:,tz+1);
    % Measure horizontal phase shifts to within +/- pi of the expected
    % angles:   
    w = w28(d,1); % expected horiz phase shift for band d
    hdphi = angle(hcp * exp(-j*w)) + w;
    %%%%%%%%% Estimation of horiz. phase gradients finishes here. %%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each subband, measure vertical phase gradients by taking the angle of summed
    % conjugate products across vertical pairs.
    prevcp = prev_complex(ty+1,:,:) .* conj(prev_complex(ty,:,:));
    refcp = ref_complex(ty+1,:,:) .* conj(ref_complex(ty,:,:));
    % The next line calculates the vertical phase gradients in subband d.
    vcp = refcp(:,tx,tz) + refcp(:,tx+1,tz) + refcp(:,tx,tz+1) + refcp(:,tx+1,tz+1)...
        + prevcp(:,tx,tz) + prevcp(:,tx+1,tz) + prevcp(:,tx,tz+1) + prevcp(:,tx+1,tz+1);
    w = w28(d,2);
    vdphi = angle(vcp * exp(-j*w)) + w;
    %%%%%%%%% Estimation of vertical phase gradients finishes here. %%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each subband, measure phase gradients between data slices by taking the angle of summed
    % conjugate products of adjacent slices.
    prevcp = prev_complex(:,:,tz+1) .* conj(prev_complex(:,:,tz));
    refcp = ref_complex(:,:,tz+1) .* conj(ref_complex(:,:,tz));
    % The next line calculates the adjacent phase gradients in subband d.
    zcp = refcp(ty,tx,:) + refcp(ty+1,tx,:) + refcp(ty,tx+1,:) + refcp(ty+1,tx+1,:) ...
        + prevcp(ty,tx,:) + prevcp(ty+1,tx,:) + prevcp(ty,tx+1,:) + prevcp(ty+1,tx+1,:); %zcp means complex in z direction
    w = w28(d,3);
    zdphi = angle(zcp * exp(-j*w)) + w;
    %%%% Estimation of phase gradients between adjacent slices finishes here. %%%%%%%

    % Measure temporal phase differences between the ref subband and prev
    % subband, by taking the angle of summed conjugate products over 2*2*2
    % regions.
    tcp = prev_complex .* conj(ref_complex);
    tcp = tcp(ty,tx,tz) + tcp(ty+1,tx,tz) + tcp(ty,tx+1,tz) + tcp(ty+1,tx+1,tz) ...
           + tcp(ty,tx,tz+1) + tcp(ty+1,tx,tz+1) + tcp(ty,tx+1,tz+1) + tcp(ty+1,tx+1,tz+1);
    tdphi = angle(tcp); % Expected phase shifts are zero in this case.
    
    % Calc. scale factors Cmag according to magnitudes of tcp, ref_complex
    % and prev_complex
    magsq = abs(ref_complex).^3 + abs(prev_complex).^3;
%    Cmag_error_term = 0.1;
%    Cmag_denominator = Cmag_error_term + magsq(ty,tx,tz) + magsq(ty+1,tx,tz) + magsq(ty,tx+1,tz) + magsq(ty+1,tx+1,tz) + magsq(ty,tx,tz+1) + magsq(ty+1,tx,tz+1) + magsq(ty,tx+1,tz+1) + magsq(ty+1,tx+1,tz+1);
    Cmag_denominator = magsq(ty,tx,tz) + magsq(ty+1,tx,tz) + magsq(ty,tx+1,tz) + magsq(ty+1,tx+1,tz) + magsq(ty,tx,tz+1) + magsq(ty+1,tx,tz+1) + magsq(ty,tx+1,tz+1) + magsq(ty+1,tx+1,tz+1);
    
    Cmag = abs(tcp).^2 ./ (Cmag_denominator + 10^(-3)); % 10^(-3) is added to the denominator, such that denominator cannot be zero
    
    % Generate constraint vectors, with units of spatial shift normalised
    % to half the image size
    % note that shift goes from -1 to +1 across the image
    cvecs{d} = cat(4,vdphi*(s2/2).*Cmag,hdphi*(s1/2).*Cmag,zdphi*(s3/4).*Cmag,tdphi.*Cmag); % Recall cvecs{d} = (y,x,z,cvec element), where d is the subband number
    if mode(1) >= it % If the first few iterations use global affine model, then assign the outermost ring of cvecs to be zeros.
            cvecs{d}([1 end],:,:,:) = 0;
            cvecs{d}(:,[1 end],:,:) = 0;
            cvecs{d}(:,:,[1 end],:) = 0;
    end
    if mode(2) == 1
       cvecs{d}([1 end],:,:,:) = 0;
       cvecs{d}(:,[1 end],:,:) = 0;
       cvecs{d}(:,:,[1 end],:) = 0;
    end
end 

%    Remove the outer parameters of constraints to reduce edge effects.
    if mode(2) == 2
        % Remove the bands nearly parallel to y axis
        cvecs{1}(:,[1 end],:,:) = 0;
        cvecs{2}(:,[1 end],:,:) = 0;
        cvecs{3}(:,[1 end],:,:) = 0;
        cvecs{4}(:,[1 end],:,:) = 0;
        cvecs{13}(:,[1 end],:,:) = 0;
        cvecs{14}(:,[1 end],:,:) = 0;
        cvecs{15}(:,[1 end],:,:) = 0;
        cvecs{16}(:,[1 end],:,:) = 0;
        cvecs{17}(:,[1 end],:,:) = 0;
        cvecs{18}(:,[1 end],:,:) = 0;
        cvecs{19}(:,[1 end],:,:) = 0;
        cvecs{20}(:,[1 end],:,:) = 0;
        
        % Remove the bands nearly parallel to x axis
        cvecs{5}([1 end],:,:,:) = 0;
        cvecs{6}([1 end],:,:,:) = 0;
        cvecs{7}([1 end],:,:,:) = 0;
        cvecs{8}([1 end],:,:,:) = 0;
        cvecs{13}([1 end],:,:,:) = 0;
        cvecs{14}([1 end],:,:,:) = 0;
        cvecs{15}([1 end],:,:,:) = 0;
        cvecs{16}([1 end],:,:,:) = 0;
        cvecs{21}([1 end],:,:,:) = 0;
        cvecs{22}([1 end],:,:,:) = 0;
        cvecs{23}([1 end],:,:,:) = 0;
        cvecs{24}([1 end],:,:,:) = 0;        
        
        % Remove the bands nearly parallel to z axis
        cvecs{1}(:,:,[1 end],:) = 0;
        cvecs{2}(:,:,[1 end],:) = 0;
        cvecs{3}(:,:,[1 end],:) = 0;
        cvecs{4}(:,:,[1 end],:) = 0;
        cvecs{5}(:,:,[1 end],:) = 0;
        cvecs{6}(:,:,[1 end],:) = 0;
        cvecs{7}(:,:,[1 end],:) = 0;
        cvecs{8}(:,:,[1 end],:) = 0;
        cvecs{9}(:,:,[1 end],:) = 0;
        cvecs{10}(:,:,[1 end],:) = 0;
        cvecs{11}(:,:,[1 end],:) = 0;
        cvecs{12}(:,:,[1 end],:) = 0;
    end
    
return;