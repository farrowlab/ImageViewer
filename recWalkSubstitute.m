function [IDX, BLOCKED, newBLOCKED] = recWalkSubstitute(M, BLOCKED, x, y, xc_theshold, maxSize)
%%
% B     =>  Basis Map size of MIP of M
% M     =>  Original Movie
% IDX   =>  2 point indexes to central pixel
    % Check Arguments
    if nargin < 6
        maxSize = length(M(:));
    end

    % set first field visited
    newBLOCKED = zeros(size(BLOCKED,1),size(BLOCKED,2));
    IDX = [x y];
    BLOCKED(x,y) = 1; newBLOCKED(x,y) = 1;
    
    GROUP = squeeze(M(IDX(1,1), IDX(1,2),:));
    
    % add 4 neighbors in todo
    TODO = [x-1 y
            x y-1
            x+1 y
            x y+1];
    while ~isempty(TODO) && sum(newBLOCKED(:)) < maxSize
        
        NEW_TODO = [];
        for i=1:size(TODO,1)
            % Check every pixel in TODO
            x = TODO(i,1);
            y = TODO(i,2);
            % Check if pixel is out of bounds
            if x <= 0 || x > size(M,1)
                continue
            end
            if y <= 0 || y > size(M,2)
                continue
            end
            % Check if we visited this pixel already
            if BLOCKED(x,y) == 1
                continue
            end
            % Block this pixel as visited
            BLOCKED(x,y) = 1; newBLOCKED(x,y) = 1;
            pixelVal = squeeze(M(x, y,:));
            
            % Test Correlation: Future add Mutual Information
            xc = corr(GROUP/size(IDX,1), pixelVal,'type','Pearson');     % 'Kendall' 'Spearman' Pearson normalize group signal     
            
            % Check if this pixel should be added to our group
            if xc > xc_theshold
                % yes, accept this pixel
                IDX = [IDX ; [x y]];
                GROUP = GROUP + pixelVal; % ????????????? does that mean that at the next iteration, the correlation test will be done on the WHOLE object, and not on the neighbouring pixel ?
                % add the neighbors of this pixel to new todolist
                NEW_TODO = [NEW_TODO;  [x-1 y
                                        x y-1
                                        x+1 y
                                        x y+1]];
            end
        end
        TODO = NEW_TODO;
    end 