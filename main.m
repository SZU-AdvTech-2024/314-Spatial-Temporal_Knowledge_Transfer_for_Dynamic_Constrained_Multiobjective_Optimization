clc
clear
close all
warning('off')
con=configure();
repeatMax=con.repeat;
functions=con.TestFunctions;
T_parameter=con.T_parameter;
popSize=con.popSize;

for rep=1:1%repeatMax
    for testFuncNo=1:size(functions,2)
        Problem=eval(functions{testFuncNo});
        for group=1:size(T_parameter,1)
            MaxIt=T_parameter(group,2);
            fprintf('\n DCF %d --- dec:%d runing on: %s, configure: %d, environment:',testFuncNo,con.dec,functions{testFuncNo},group);
            res=DCMFEA(Problem,popSize,MaxIt,T_parameter,group);
            [igd]=computeMetrics(res);
            fprintf('\n igd:%d',mean(igd));
        end 
    end
end

function [IGD_T]=computeMetrics(resStruct)
for T=1:size(resStruct,2)
    cpof=resStruct{T}.CPOF;
    cpof(imag(cpof)~=0) = abs(cpof(imag(cpof)~=0));
    pof=resStruct{T}.POF;
    pof(imag(pof)~=0) = abs(pof(imag(pof)~=0));
    if size(cpof,2)==0
        IGD_T(T)=IGD(pof',resStruct{T}.turePOF);
    else
        IGD_T(T)=IGD(cpof',resStruct{T}.turePOF);
    end
end
resIGD=mean(IGD_T);
end
