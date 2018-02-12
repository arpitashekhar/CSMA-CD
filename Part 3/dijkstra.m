
function [costs, paths] = dijkstra(Edges, Weight)
    % Process inputs
    n = size(Edges, 1);
    source_ids = (1:n);
    dest_ids = (1:n);

    [E, cost] = adj2edge(Edges, Weight);
    E = E(:,[2 1]);
       
    % Initialize output variables
    L = length(source_ids);
    M = length(dest_ids);
    costs = zeros(L,M);
    paths = num2cell(NaN(L,M));
    
    % Find the minimum costs and paths using Dijkstra's Algorithm
    for k = 1:L
        
        % Initializations
        iTable = NaN(n,1);
        minCost = Inf(n,1);             % initializing minimum cost values for every destination to infinity
        isMin = false(n,1);             % flag for identifying if we have found the shortest path to that node
        path = num2cell(NaN(n,1));
        I = source_ids(k);              % node at index k
        minCost(I) = 0;
        iTable(I) = 0;
        isMin(I) = true;
        path(I) = {I};
        
        % Execute Dijkstra's Algorithm for this vertex
        while any(~isMin(dest_ids))
            
            % Update the table
            jTable = iTable;
            iTable(I) = NaN;
            nodeIndex = find(E(:,1) == I);          % finding the adjacent nodes
            
            % Calculate the costs to the adjacent nodes and record paths
            for x = 1:length(nodeIndex)
                J = E(nodeIndex(x),2);
                if ~isMin(J)
                    c = cost(I,J);
                    empty = isnan(jTable(J));
                    if empty || (jTable(J) > (jTable(I) + c))      % check if the new path has a lower cost than before
                        iTable(J) = jTable(I) + c;
                        path{J} = [path{I} J];
                        
                    else
                        iTable(J) = jTable(J);
                    end
                end
            end
            
            % Find values in the table
            K = find(~isnan(iTable));
            if isempty(K)
                break
            else
                % Set the minimum value for the nodes to true
                [~,N] = min(iTable(K));
                I = K(N);
                minCost(I) = iTable(I);
                isMin(I) = true;
            end
        end
        
        % Store costs and paths
        costs(k,:) = minCost(dest_ids);
        paths(k,:) = path(dest_ids);
    end          
end

% Convert adjacency matrix to edge list
function [E,C] = adj2edge(Edges, Weight)
    [I,J] = find(Edges);
    E = [I J];
    C = Weight;
end