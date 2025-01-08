function res=DCMFEA(Problem,popSize,MaxIt,T_parameter,group)
Nfor2=20;
    F(1) = struct('cdata',[],'colormap',[]);
    rmp=1;
    muc = 1; 
    mum = 1;
    prob_vswap = 0; 
    for T=1:T_parameter(group,3)/T_parameter(group,2)
        fprintf(' %d',T);
        t = 2/T_parameter(group,1)*(T-1);  
        mainTask=Problem;
        auxiliaryTask=creatAuxiliaryTask(Problem);
        if T==1
            [POF,population]=MOMFEA_v1(t,mainTask,popSize-Nfor2,auxiliaryTask,Nfor2,rmp,MaxIt,muc,mum,prob_vswap);  
        else

          [POF,population]=MOMFEA_v1(t,mainTask,popSize-Nfor2,auxiliaryTask,Nfor2,rmp,MaxIt,muc,mum,prob_vswap,population);  
        end
        
        
        CPOF=getFesiableSolutions(t,POF,mainTask);

        res{T}.turePOF=Problem.getCPF(400,t);  
        res{T}.POF=POF;
        res{T}.CPOF=CPOF;
    end
end




function testInstance=creatAuxiliaryTask(Problem)
    testInstance.getFesiableRegion=Problem.getFesiableRegion;
    testInstance.getUPF=Problem.getUPF;
    testInstance.getCPF=Problem.getCPF;
    testInstance.getObj=@getObj;
    testInstance.getBound=Problem.getBound;
    testInstance.M=2;
    testInstance.D=10;
    testInstance.LBound=zeros(1,testInstance.D);
    testInstance.UBound=ones(1,testInstance.D);

end

function PopObj=getObj(X,M,D,t,Problem)
    PopObj1=Problem.getObj(X,M,D,t);
    PopObj2=Problem.getCV(PopObj1,t);
    PopObj=[PopObj1];
end

