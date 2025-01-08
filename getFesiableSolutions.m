
function CPF=getFesiableSolutions(t,POF,Problem)
    CV=Problem.getCV(POF',t);
    for i=1:size(CV,2)
        index= CV(:,i)>0;
        POF(:,index)=[];
        CV(index,:)=[];
    end
    CPF=POF;
end