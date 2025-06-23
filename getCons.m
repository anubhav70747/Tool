function Cons = getCons(GammaRand,Gamma)
    tolerance = 10^-7;
    [indicesList,zeroIndecesList] = findZerosAndSimilarEntries(GammaRand, tolerance);
    numel(indicesList)
    Cons = [];
    zeroIndecesList
    for i =1:length(zeroIndecesList)
        zeroIndecesList{i}
        [i,j] = zeroIndecesList{i};
        Cons = [Cons; Gamma(i,j)==0];
    end
    length(indicesList)
    for i =1:length(indicesList)
        listSimilarIndices=indicesList{i};
        [r,c] = size(listSimilarIndices);

        if r>1
            
            for j = 2:r
                Cons=[Cons;Gamma(listSimilarIndices(1,1),listSimilarIndices(1,2))==Gamma(listSimilarIndices(j,1),listSimilarIndices(j,2))];
            end
        end
    end
    
    
    
