function adjMatrix = RandomTree(numberNode,Flag)
adjMatrix = zeros(numberNode);
%% Random tree
if Flag == 0
    for ii = 2:numberNode
        k = ii-1;
        adjMatrix(ii,k) = 1;
        adjMatrix(k,ii) = 1;
    end
else
    for ii = 2:numberNode
        k = 1;
        adjMatrix(ii,k) = 1;
        adjMatrix(k,ii) = 1;
    end
end
