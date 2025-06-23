function [S,PXYAB]=BellTool(X,Y,A,B)
    gammaR=RandomMomentMatrixLevel3(X,Y,A,B);
    l=length(gammaR);
    gamma = sdpvar(l,l,'hermitian','real');
    con=[];
    con = [con;gamma(1,1)==1];
    con = [con;gamma>=0];
    con=[con;getCons(gammaR,gamma)];
    idxA = @(x,a) (A-1)*(x-1)+a+1;
    idxB = @(y,b) (A-1)*X+1 + (B-1)*(y-1)+b;
    PXA=cell(X,A);
    for x=1:X
        sum=0;
        for a=1:A-1
            PXA{x,a}=gamma(1,idxA(x,a));
            sum=sum+PXA{x,a};
        end
        PXA{x,A}=1-sum;
       % con=[con;PXA{x,A}>=0]
    end

    PYB=cell(Y,B);
    for y=1:Y
        sum=0;
        for b=1:B-1
            PYB{y,b}=gamma(1,idxB(y,b));
            sum=sum+PYB{y,b};
        end
        PYB{y,B}=1-sum;
       % con=[con;PYB{y,B}>=0]
    end
    PXYAB=cell(X,Y,A,B);
    for x=1:X
        for y= 1:Y
            for a=1:A-1
                for b=1:B-1
                    PXYAB{x,y,a,b}=gamma(idxA(x,a),idxB(y,b));
                end
            end
        end
    end


    S=[];
S1=real(-PYB{1,1}-2*PXA{1,1}-PXA{2,1}+PXYAB{1,1,1,1}+PXYAB{2,1,1,1}+PXYAB{3,1,1,1}+PXYAB{1,2,1,1}+PXYAB{2,2,1,1}-PXYAB{3,2,1,1}+PXYAB{1,3,1,1}-PXYAB{2,3,1,1});
   
   diagnostics = optimize([con], -S1, sdpsettings('solver', 'mosek'));
    S = [S;value(S1)];