function av = estshhemm33( w, nit, levelsel, ...
     avlevel, mode, qscale, sizeqfilt ) 
% ESTSHHEMM33 
% This version of estshhemm generates av (affine parameters) at the 
% resolution of 'avlevel'. av will be used to generate full resolution
% motion field later to register the datasets.
% xa is the dataset to be registered and xb is the reference dataset.
% w is the expected phase rotations per sample.
% nit is the number of iterations performed
% levelsel specifies the range of the CWT levels used at each iteration.
% avlevel is the level of the DTCWT subbands which have the same resolution
% as the avec matrix
% sizeqfilt is the size of the Q matrix smoothing filter (no. of taps = 2*sizeqfilt +1).
% shift_x_tot, shift_y_tot and shift_y_tot are 3-D arrays which describe
% the motion vector in x, y and z direction respectively.
% delta(it) gives the mean error of at each iteration.
% av is the affine vectors which will then be used to generate the shift
% vectors for image registration.

global ya ya_before_reg % ya and ya_before_reg are the shifted highpass bands of xa and original highpass bands of xa respectively
global band_size
if nargin < 7, sizeqfilt = 2; end %sizeqfilt is the size of the Q matrix smoothing filter
if nargin < 6, qscale = 1; end  %The scale factor for the quiver plots
if nargin < 5; mode = 0; end
%if length(mode) < 5, mode = [mode 0 0 0 0]; end


fig = gcf +1;

% Check if the matrices of the LLL_band and highpass bands at each level to
% be of the same size
% if (size(xa) ~= size(xb))
%     error ('Input matrices xa and xb do not have the same size');
% end

% Def. the expected phase shifts along dim1 dim2 and dim3 for each of the 28
% subbands.
% w28 is a matrix with 28 rows and 3 columns, each row specifies the
% expected phase shifts along 3 dimensions (x, y and z)
% w28 = [HLL; LHL; HHL; LLH; HLH; LHH; HHH];
w28 =  [w(2) w(1) w(1); -w(2) w(1) w(1); w(2) -w(1) w(1); -w(2) -w(1) w(1);  %HLL sub-band
        w(1) w(2) w(1); -w(1) w(2) w(1); w(1) -w(2) w(1); -w(1) -w(2) w(1);  %LHL sub-band
        w(2) w(2) w(1); -w(2) w(2) w(1); w(2) -w(2) w(1); -w(2) -w(2) w(1);  %HHL sub-band
        w(1) w(1) w(2); -w(1) w(1) w(2); w(1) -w(1) w(2); -w(1) -w(1) w(2);  %LLH sub-band
        w(2) w(1) w(2); -w(2) w(1) w(2); w(2) -w(1) w(2); -w(2) -w(1) w(2);  %HLH sub-band
        w(1) w(2) w(2); -w(1) w(2) w(2); w(1) -w(2) w(2); -w(1) -w(2) w(2);  %LHH sub-band
        w(2) w(2) w(2); -w(2) w(2) w(2); w(2) -w(2) w(2); -w(2) -w(2) w(2)]; %HHH sub-band
   
 %Def length of the avec = 12
 lenavec = 12;
 
 srefh = band_size(avlevel,:); %srefh contains the size of the sub-band at level = avlevel
 s2 = srefh(1)/2; 
 s1 = srefh(2)/2;
 s3 = srefh(3); % Now s2,s1,s3 give the size of 1 of the 28 subbands
 
 % Initialise affine vector 
 av = zeros(s2-1, s1-1, s3/2-1, lenavec); % affine vector is of the size of 1 of the 28 subband minus 1 in each direction (except in z direction, frames are scaled down by a factor of 2)

 
 % Perform nit iterations on shift and del so that del is more accurate
 for it = 1:nit,
      levels = levelsel(min(it,size(levelsel,1)),:);
      qvecsum = 0;
      ya = ya_before_reg;
      % Loop for each level at the current iteration.
      for lev = levels(1):-1:levels(2)

          % Calc. size of subband at lev
          sref_subband = band_size(lev,:); %sref_subband is the size of 1 of the 8 subbands at current level
          s2 = sref_subband(1)/2;
          s1 = sref_subband(2)/2;
          s3 = sref_subband(3);
           
          [shift_x, shift_y, shift_z] = affineshift3(av,lev,avlevel); %The obtained shift vectors has the same size as that of the subband (except the frame no is reduced by a factor of 2)
          % Use shift (adjusted to voxel units instead of half-image units) to
          % shift xa.
          method = 'cub';
          % xa = xa_before_reg;
          % ya = ya_before_reg;
          
          yi_distance = 2/band_size(1,1)*(2^lev);
          xi_distance = 2/band_size(1,2)*(2^lev);
          zi_distance = 2/band_size(1,3)*(2^lev);
          % Use shift vectors (adjusted for pel units instead of half-image units) to shift subband.
          shift_xa_cwt_bands3('hi', lev,shift_x/xi_distance,shift_y/yi_distance,shift_z/zi_distance,method); % This line shifts the 28 subbands at the specified level
  
          % Calc. motion constraints at each point. 
          % cvecs is a 5-D array in (y,x,z,band,cvec element), where the
          % number of bands is 28
          cvecs = zeros(s2-1,s1-1,s3/2-1,28,lenavec+1); %Initialise cvecs, noting that although cvec has 4 elements for a location, cvecs will be increased later by multiplying K_hat
          
          cvecs_temp = constraint3(lev, w28, mode, it); %temporary cvecs which contains only 4 elements for a location
          
          for d=1:28
              cvecs(:,:,:,d,[1 2 3 lenavec+1]) = cvecs_temp{d};
          end
          
          % The next section multiplies cvecs_temp by the K_hat matrix transpose
          % K matrix is [1 0 0 y 0 0 x 0 0 z 0 0; 0 1 0 0 y 0 0 x 0 0 z 0; 0 0 1 0 0 y 0 0 x 0 0 z]
          % K_hat matrix is [K,[0 0 0]';[0 0 0 0 0 0 0 0 0 0 0 0 1]];

          % yi, xi and zi go from -1 to +1 across the subband.
          for k = 1:(s2-1)
              yi = (k-s2/2)*yi_distance;
              cvecs(k,:,:,:,[4 5 6]) = cvecs(k,:,:,:,[1 2 3]) * yi;
          end
          
          for k = 1:(s1-1)
              xi = (k-s1/2)*xi_distance;
              cvecs(:,k,:,:,[7 8 9]) = cvecs(:,k,:,:,[1 2 3]) * xi;
          end
          
          for k = 1:(s3/2-1)
              zi = (k-s3/4)*zi_distance;
              cvecs(:,:,k,:,[10 11 12]) = cvecs(:,:,k,:,[1 2 3]) * zi;
          end
          % Now cvecs contains the term K_hat' * c
          
          % Form the Q-matrices as vectors of the lower left half of the
          % symmetric Q-matrix. Q = (K_hat'*c*c'*K_hat)
          scv = size(cvecs); 
          c5 = scv(5); % this c5 also equals to lenavec+1 = 13
          qvecs = zeros([scv(1:3) c5*(c5+1)/2]);
          iQ = zeros(c5,c5);
          k = 1;
          for k1 = 1:c5
              for k2 = k1:c5
                  qvecs(:,:,:,k) = sum(cvecs(:,:,:,:,k1) .* cvecs(:,:,:,:,k2),4); % Sum over all 28 subbands.
                  iQ(k1,k2) = k; iQ(k2,k1) = k;
                  k=k+1;
              end
          end       

          % Change the sample rate of qvecs, back to that of av,
          % remembering that the sizes are always odd and one less than the
          % subband sizes.          
          if lev < avlevel
              %Reduce the sample rate.
              for i = 1:(avlevel-lev)
                  % check if the coarser level subband is exactly twice as big as the
                  % finer level subband. If not, do appending.
                  if band_size(lev+i-1,1)/4 ~= band_size(lev+i,1)/2 
                     qvecs = cat(1,qvecs(1,:,:,:),qvecs,qvecs(end,:,:,:));
                  end
                  if band_size(lev+i-1,2)/4 ~= band_size(lev+i,2)/2 
                     qvecs = cat(2,qvecs(:,1,:,:),qvecs,qvecs(:,end,:,:));
                  end
                  if band_size(lev+i-1,3)/4 ~= band_size(lev+i,3)/2 
                     qvecs = cat(3,qvecs(:,:,1,:),qvecs,qvecs(:,:,end,:));
                  end
                  

                  % First reduce no of rows, using [0.25 0.5 0.25] filter.
                  t1 = 2:2:(size(qvecs,1)-1);
                  qvecs = 0.25*qvecs(t1-1,:,:,:) + 0.5*qvecs(t1,:,:,:) + 0.25*qvecs(t1+1,:,:,:);
                  % Now reduce no of columns.
                  t2 = 2:2:(size(qvecs,2)-1);
                  qvecs = 0.25*qvecs(:,t2-1,:,:) + 0.5*qvecs(:,t2,:,:) + 0.25*qvecs(:,t2+1,:,:);
                  % Now reduce no of frames.
                  t3 = 2:2:(size(qvecs,3)-1);
                  qvecs = 0.25*qvecs(:,:,t3-1,:) + 0.5*qvecs(:,:,t3,:) + 0.25*qvecs(:,:,t3+1,:);
                  
              end
          elseif lev > avlevel,
              for i = 1:(lev-avlevel)
                  % First increase no of rows, using [0.5 1 0.5] filter.
                  sr = size(qvecs,1);
                  qv = [];
                  qv([0:sr]*2+1,:,:,:) = 0.5 * (qvecs([1 1:sr],:,:,:) + qvecs([1:sr sr],:,:,:)); 
                  qv([1:sr]*2,:,:,:) = qvecs; 
                  qvecs = qv;
                  % Now increase no of columns
                  sc = size(qvecs,2);
                  qv = [];
                  qv(:,[0:sc]*2+1,:,:) = 0.5 * (qvecs(:,[1 1:sc],:,:) + qvecs(:,[1:sc sc],:,:)); 
                  qv(:,[1:sc]*2,:,:) = qvecs; 
                  qvecs = qv;
                  % Now increase no of frames
                  sf = size(qvecs,3);
                  qv = [];
                  qv(:,:,[0:sf]*2+1,:) = 0.5 * (qvecs(:,:,[1 1:sf],:) + qvecs(:,:,[1:sf sf],:)); 
                  qv(:,:,[1:sf]*2,:) = qvecs; 
                  qvecs = qv;
                  
                  % Check if the coarser level subband is exactly twice as big as
                  % the finer level subband. If not, discard 2 outer most
                  % coefs.
                  if band_size(lev-i,1)/4 ~= band_size(lev-i+1,1)/2
                      qvecs = qvecs(2:end-1,:,:,:);
                  end
                  if band_size(lev-i,2)/4 ~= band_size(lev-i+1,2)/2
                      qvecs = qvecs(:,2:end-1,:,:);
                  end
                  if band_size(lev-i,3)/4 ~= band_size(lev-i+1,3)/2
                      qvecs = qvecs(:,:,2:end-1,:);
                  end         
                  
              end
          end
          qvecsum = qvecsum + qvecs; % Accumulate qvecs across levels of CWT.
      end % End of loop for each level of CWT
       
    qtot = squeeze(sum(sum(sum(qvecsum))))/size(qvecsum,1)/size(qvecsum,2)/size(qvecsum,3);
    normalised_qtot = qtot/norm(qtot);
    
    
    if mode(1) < it
      % Filter the q vectors with a linear ramp filter
      shq = sizeqfilt(min(it,end)) + 1; %  % N.B sizeqfilt is the size of the Q matrix smoothing filter
      hq = [1:shq  (shq-1):-1:1].'; 
      hq = hq/sum(hq);
       for k=1:size(qvecsum,4) %i.e c5*(c5+1)/2
              qvecsum_k = qvecsum(:,:,:,k);

              sqvecsum_k = size(qvecsum_k);
              s1 = 1:sqvecsum_k(1); s2 = 1:sqvecsum_k(2); s3 = 1:sqvecsum_k(3);
              for f = s2,
                  y = reshape(qvecsum_k(s1,f,s3),sqvecsum_k([1 3])).';
                  % Do filtering on frames
                  qvecsum_k(s1,f,s3) = colfilter(y,hq).';
              end

              for f = s3,
                  % Do filtering on rows and then columns
                  qvecsum_k(s1,s2,f) = colfilter(colfilter(qvecsum_k(s1,s2,f),hq).', hq).';
              end
              
              qvecsum(:,:,:,k) = qvecsum_k;
       end
    elseif mode(1) >=it
        for k=1:size(qvecsum,4)
            qvecsum(:,:,:,k) = qtot(k);
        end
    end



      % Now extract Q matrix, q vector and q0 scalar from each qvec vector, and
      % solve for the affine vector, a, which minimises the constraint
      % error energy.
      iQmat = iQ(1:(end-1),1:(end-1));  % Indices to convert qvec to Q.
      iqv = iQ(1:(end-1),end);  % Indices to convert qvec to q.
      iq0 = iQ(end,end);  % Index to convert qvec to q0.
      sqv = size(qvecsum);
      del = zeros(sqv(1:3));
      avec = zeros([sqv(1:3) lenavec]); % initialise the avec, avec is a 4-D array [y,x,z,avec elements]
      qf = zeros(sqv(4),1);
      

      % The information propagation scheme:
      % To use the information propagation scheme, uncomment the next line,
      % and assign the lamda in the next section of for loop to be zero.
      % qvecsum = info_propagation(qvecsum);
      
      % Loops for each Q matrix locality.
      for column = 1:sqv(2)
          for row = 1:sqv(1)
              for frame = 1:sqv(3)
                  qf(:) = qvecsum(row,column,frame,:);
                  
                  lamda = 0.1*(norm(qf)+1); %1 is added to avoid the case when norm(qf) is zero
                  %lamda = 0;
                  % Now modify qf to be the local qf plus a small portional
                  % of the global qf
                  qf = qf + lamda*normalised_qtot;
                  
                  % Extract Qmat, qv and q0.
                  Qmat = qf(iQmat);
                  qv = qf(iqv);
                  q0 = qf(iq0);
                  % Solve for affine parameters: 
                  a = -pinv(Qmat) * qv;
                  avec(row,column,frame,:) = a;
                  % Calculate the error energy:
                  del(row,column,frame) = a.'*Qmat*a + 2*qv.'*a + q0;
              end
          end
      end
      
      % Update av.
      av = av + avec;
    
      % Create shift vector at same resolution as av.
      sav = size(av) + [1 1 1 0];
      % Obtain the extra shift in x, y and z directions on current
      % iteration
      [shift_x_current_it, shift_y_current_it, shift_z_current_it] = affineshift3(avec, avlevel, avlevel);
      yi_distance = 2/band_size(1,1)*(2^avlevel);
      xi_distance = 2/band_size(1,2)*(2^avlevel);
      zi_distance = 2/band_size(1,3)*(2^avlevel);
      shift_mag = sqrt((shift_x_current_it/xi_distance).^2 + (shift_y_current_it/yi_distance).^2 + (shift_z_current_it/zi_distance).^2)*(2^avlevel);
      delta = mean(del(:));
      max_shift_it = max(max(max(shift_mag)));
	  mean_shift_current_it = sum(sum(sum(shift_mag)))/size(shift_mag,1)/size(shift_mag,2)/size(shift_mag,3);  % Display extra shift magnitude for this iteration. 
      % Monitor the mean squared error and the mean affine vector (per
      % sample).
      string = sprintf('Iteration %d',it);
      disp(string)
      string = sprintf('Operating on levels from %d to %d',levels(1),levels(2));
      disp(string)
      string = sprintf('Squared error of the model: %f',delta);
      disp(string)
      string = sprintf('The max of extra shifts in the current iteration: %f', max_shift_it);
      disp(string)
      string = sprintf('The average of extra shifts in the current iteration: %f \n', mean_shift_current_it);
      disp(string)
      
     % Update shift.
     [shift_x_tot shift_y_tot shift_z_tot] = affineshift3(av,avlevel,avlevel);
      %{
      % Open figure to plot 3-D shift vectors
      figure(fig); set(gcf, 'DefaultLineLineWidth',1.5);
      % Generate the x, y and z points at which the shift vectors are to be
      % plotted
      z_axis = zeros(sav(1:3));
      for k=1:sav(3)
          for t=1:sav(1)
              x_axis(t,:,k) = [1:sav(2)];
          end
          for t=1:sav(2)
              y_axis(:,t,k) = [1:sav(1)]';
          end
          z_axis(:,:,k) = k;
      end
      % Now plot the shift vectors
      quiver3(x_axis,y_axis,z_axis,-shift_x_tot/xi_distance*qscale, -shift_y_tot/yi_distance*qscale, -shift_z_tot/zi_distance*qscale,0);
      set(gcf,'position',[700 20 500 672]);
      set(gca,'position',[0.01 0.01 .98 .98]);
      axis equal;
      axis ij;
      axis([-1 sav(2)+2 -1 sav(1)+2 -1 sav(3)+2]); % set x,y,z axes      
      settitle(sprintf('%d iteration, output shift vectors',it));
      drawnow
     %}
      fig = fig + 1;

      
 end % End of current iteration
 
 return;