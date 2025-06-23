function Gamma = RandomMomentMatrixLevel3(X, Y, A, B)
% RANDOMMOMENTMATRIXLEVEL3  Build the level-3 NPA moment matrix Γ.
%
%   Γ = RandomMomentMatrixLevel3(X, Y, A, B) returns the real
%   moment matrix for a bipartite Bell‐type scenario in which
%     • Alice has X measurement settings with A outcomes each;
%     • Bob   has Y measurement settings with B outcomes each.
%
%   Projective measurements and the global state are sampled at random
%   (using the helper routines RandomProjectiveMeasurements and
%   RandomStateVector that you already have).  The operator list contains
%   every monomial of length ≤ 3 in Alice/Bob projectors (with the last
%   outcome of every setting omitted because of the “sum-to-identity”
%   relation).
%
%   The code is written for MATLAB / Octave R2020b+ and avoids anonymous
%   functions that assign to outer-scope variables, using nested helper
%   functions instead.

% ------------------------------------------------------------------------
% 1.  random POVMs and random pure state
% ------------------------------------------------------------------------
[AliceM, da] = RandomProjectiveMeasurements(X, A,2);   % cell{X,A}
[BobM,   db] = RandomProjectiveMeasurements(Y, B,2);   % cell{Y,B}

psi = RandomStateVector(da * db);
rho = psi * psi';                                      % |ψ⟩⟨ψ|

% ------------------------------------------------------------------------
% 2.  build the operator list (monomials of length ≤ 3)
% ------------------------------------------------------------------------
OperatorList = { eye(da * db) };    % length-0 (identity)

% -- helper functions that append to OperatorList (must be nested) -------
    function addA(Aop)
        OperatorList{end+1} = kron(Aop, eye(db));
    end
    function addB(Bop)
        OperatorList{end+1} = kron(eye(da), Bop);
    end
    function addAB(Aop, Bop)
        OperatorList{end+1} = kron(Aop, Bop);
    end
% ------------------------------------------------------------------------

% ----- length-1 monomials ----------------------------------------------
for x = 1:X
    for a = 1:A-1                    % omit last projector
        addA(AliceM{x, a});
    end
end
for y = 1:Y
    for b = 1:B-1
        addB(BobM{y, b});
    end
end

% ----- length-2 monomials ----------------------------------------------
% AA pairs
for x1 = 1:X
    length(OperatorList)
    for a1 = 1:A-1
        for x2 = 1:X
            for a2 = 1:A-1
                if x1 ~= x2
                  addA(AliceM{x1, a1} * AliceM{x2, a2});
                end
            end
        end
            for y = 1:Y
            for b = 1:B-1
                addAB(AliceM{x1,a1},BobM{y,b})
            end
        
        end
    end
end
% AB pairs
 for x = 1:X
    for a = 1:A-1
        for y = 1:Y
            for b = 1:B-1
            %    addAB(AliceM{x,a},BobM{y,b})
            end
        end
    end
 end

 % BB pairs
for y1 = 1:Y
    for b1 = 1:B-1
        for y2 = 1:Y
            for b2 = 1:B-1
                if y1 ~= y2
                    addB(BobM{y1, b1} * BobM{y2, b2});
                end
            end
        end
    end
end

% AAA
for x1 = 1:X
    for a1 = 1:A-1
        for x2 = 1:X
            for a2 = 1:A-1
                for x3 = 1:X
                    for a3 = 1:A-1
                        if x1~=x2 && x2~=x3 
                          addA(AliceM{x1,a1} * AliceM{x2,a2} * AliceM{x3,a3});
                        end
                    end
                end
                for y = 1:Y
            for b = 1:B-1
                        if x1~=x2 
                            addAB(AliceM{x1,a1} * AliceM{x2,a2},BobM{y,b})
                        end
                    end
        
                end
            end
        end
        for y1 = 1:Y
            for b1 = 1:B-1
                for y2 = 1:Y
                    for b2 = 1:B-1
                        if y1 ~= y2
                            addAB(AliceM{x1,a1},BobM{y1, b1} * BobM{y2, b2});
                        end
                    end
                end
            end
        end
    end
end

for y1 = 1:Y
    for b1 = 1:B-1
        for y2 = 1:Y
            for b2 = 1:B-1
                for y3 = 1:Y
                    for b3 = 1:B-1
                        if y1~=y2 && y2~=y3 
                           addB(BobM{y1,b1} * BobM{y2,b2} * BobM{y3,b3});
                        end
                    end
                end
                
            end
        end
        
    end
end
% AA-B
for x1 = 1:X
    for a1 = 1:A-1
        for x2 = 1:X
            for a2 = 1:A-1
                for y = 1:Y
                    for b = 1:B-1
                        if x1~=x2
                   %         addAB( AliceM{x2,a2}*AliceM{x1,a1}, BobM{y,b});
                        end
                    end
                end
            end
        end
    end
    end
% A-BB
for x = 1:X
    for a = 1:A-1
        for y1 = 1:Y
            for b1 = 1:B-1
                for y2 = 1:Y
                    for b2 = 1:B-1
                        if y1~=y2
                   %         addAB(AliceM{x,a}, BobM{y1,b1} * BobM{y2,b2});
                        end
                    end
                end
            end
        end
    end
end



% ------------------------------------------------------------------------
% 3.  assemble the moment matrix Γ
% ------------------------------------------------------------------------
L = numel(OperatorList)
Gamma = zeros(L);

for i = 1:L
    Oi = OperatorList{i};
    for j = 1:L
        Gamma(i, j) = (trace(rho*Oi'* OperatorList{j}));
    end
end


end