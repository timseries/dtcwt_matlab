function qvecsum = info_propagation( qvecsum )
% function qvecsum = info_propagation(qvecsum)
% This function takes the input qvecsum and determines which localities are
% informative. We then propagate the qvecsum values of those localities 
% which are informative to the regions with no information.

sqv = size(qvecsum);
qvec_energy = zeros(sqv(1),sqv(2),sqv(3));

for column = 1:sqv(2)
    for row = 1:sqv(1)
        for frame = 1:sqv(3)
            qvec_energy(row, column, frame) = qvecsum(row, column, frame, 1); 
        end
    end
end

% Setting the threshold as the median of the qvec_energy. The threshold is
% used to decide whether a locality contains motion information or it is
% nearly ill-conditioned and should be discarded.
threshold = median(median(median(qvec_energy)));

% info_locality_map stores the locations of informative localities 
info_locality_map = (qvec_energy>=threshold);

while ~all(all(all(info_locality_map)))
    [info_row_index, info_column_index] = find(info_locality_map);
    info_frame_index = ceil(info_column_index/sqv(2));
    info_column_index = rem(info_column_index,sqv(2));
    info_column_index(info_column_index==0) = sqv(2);

    qvecsum_temp = qvecsum;
    
    for column = 1:sqv(2)
        for row = 1:sqv(1)
            for frame = 1:sqv(3)
                if info_locality_map(row,column,frame) == 0,
                    % compute the distance between the current locality and
                    % the informative localities
                    distance = abs(info_row_index-row) + abs(info_column_index-column) + abs(info_frame_index-frame);
                    % When we have found a informative locality which is the
                    % neighbouring locality of the current locality, then
                    % we figure out the row, column and frame indices of
                    % the neighbouring locality.
                    searched_info_row_index = info_row_index(distance==1);
                    searched_info_column_index = info_column_index(distance==1);
                    searched_info_frame_index = info_frame_index(distance==1);
                    
                    if isempty(searched_info_row_index) ~= 1
                        qvecsum_temp(row,column,frame,:) = zeros(size(qvecsum_temp(row,column,frame,:)));
                        for i = 1:length(searched_info_row_index)
                            % summing up the values of all nearest informative localities
                             qvecsum_temp(row,column,frame,:) = qvecsum_temp(row,column,frame,:) + qvecsum(searched_info_row_index(i),searched_info_column_index(i),searched_info_frame_index(i),:);
                        end
                        qvecsum_temp(row,column,frame,:) = qvecsum_temp(row,column,frame,:)/length(searched_info_row_index); % Take the average
                        info_locality_map(row,column,frame) = 1; % update the info_locality_map
                    end
                end
            end
        end
    end
    qvecsum = qvecsum_temp;
end

return