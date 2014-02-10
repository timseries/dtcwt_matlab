function process_dataset( action, dim, sample )
% PROCESS_DATASET
% This function takes a dataset (defined as a global variable x), append
% zeros at its edges or trim the dataset edges.
% E.g. process_dataset('trim', 'row', [1 2]) trims the top row and bottom
% two rows
global x

if strcmp(action,'trim') == 1
    if strcmp(dim,'row') == 1
        x = x(1+sample(1):end-sample(2),:,:);
    elseif strcmp(dim,'column') == 1
        x = x(:,1+sample(1):end-sample(2),:);
    elseif strcmp(dim,'frame') == 1
        x = x(:,:,1+sample(1):end-sample(2));
    else
        error('dim must be row, column or frame');
    end
elseif strcmp(action,'extend') == 1
    sx = size(x);
    if strcmp(dim,'row') == 1
       x = cat(1,zeros(sample(1),sx(2),sx(3)),x,zeros(sample(2),sx(2),sx(3)));
    elseif strcmp(dim,'column') == 1
        x = cat(2,zeros(sx(1),sample(1),sx(3)),x,zeros(sx(1),sample(2),sx(3)));
    elseif strcmp(dim,'frame') == 1
       x = cat(3,zeros(sx(1),sx(2),sample(1)),x,zeros(sx(1),sx(2),sample(2)));
    else
        error('dim must be row, column or frame');
    end
end
            