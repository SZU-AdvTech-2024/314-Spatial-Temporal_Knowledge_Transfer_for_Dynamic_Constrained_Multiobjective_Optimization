function [POF,population]=MOMFEA_v1(t,mainTask,pop1,auxiliaryTask,pop2,rmp,gen,muc,mum,prob_vswap,initPop)
[L1,U1]=mainTask.getBound(t);
[L2,U2]=auxiliaryTask.getBound(t);


pop=pop1+pop2;
if mod(pop,2)~=0
    pop=pop+1;
    pop2=pop2+1;
end

dim1=length(L1);
dim2=length(L2);
dim=max([dim1,dim2]);


store=[];
if nargin==11
    for i=1:length(initPop)
        population(i)=Chromosome;
        population(i).rnvec=initPop(i).rnvec;
        population(i).skill_factor=initPop(i).skill_factor;
    end
    for i=length(initPop)+1:pop
        population(i)=Chromosome;
        population(i)=initialize(population(i),dim);
        if i<=pop1
            population(i).skill_factor=1;
        else
            population(i).skill_factor=2;
        end
    end
else
    for i=1:pop
        population(i)=Chromosome;
        population(i)=initialize(population(i),dim);
        if i<=pop1
            population(i).skill_factor=1;
        else
            population(i).skill_factor=2;
        end
    end
end



for i=1:pop
    population(i)=evaluate(t,population(i),L1,U1,mainTask,L2,U2,auxiliaryTask,dim1,dim2);
end

%%%%%%%%%%%%5
if nargin==11
    targetpopulation_T1=population([population.skill_factor]==1);
    targetpopulation_T2=population([population.skill_factor]==2);
    sourcrpopulation_T1=initPop([initPop.skill_factor]==1);
    sourcrpopulation_T2=initPop([initPop.skill_factor]==2);
    targetpopulation_T1=NDS_CV(targetpopulation_T1,mainTask);
    [targetpopulation_T2,~]=SolutionComparison.nondominatedsort(targetpopulation_T2,length(targetpopulation_T2),2);
    for j=1:length(sourcrpopulation_T1)
        T1xsource(j,:)=sourcrpopulation_T1(j).rnvec;
        T1ysource(j)=sourcrpopulation_T1(j).front;
        T1xtarget(j,:)=targetpopulation_T1(j).rnvec;
        T1ytarget(j)=targetpopulation_T1(j).front;
    end
    for j=1:length(sourcrpopulation_T2)
        T2xsource(j,:)=sourcrpopulation_T2(j).rnvec;
        T2ysource(j)=sourcrpopulation_T2(j).front;
        T2xtarget(j,:)=targetpopulation_T2(j).rnvec;
        T2ytarget(j)=targetpopulation_T2(j).front;
    end

    for iii=1:length(targetpopulation_T1)
        targetPopObjT1(iii,:)=targetpopulation_T1(iii).objs_T1;
    end
    for iii=1:length(sourcrpopulation_T1)
        sourcePopObjT1(iii,:)=sourcrpopulation_T1(iii).objs_T1;
    end
    for iii=1:length(targetpopulation_T2)
        targetPopObjT2(iii,:)=targetpopulation_T2(iii).objs_T2;
    end
    for iii=1:length(sourcrpopulation_T2)
        sourcePopObjT2(iii,:)=sourcrpopulation_T2(iii).objs_T2;
    end
    
    targetT1CrowdDis = CrowdingDistance(targetPopObjT1,T1ytarget);
    targetT2CrowdDis = CrowdingDistance(targetPopObjT2,T2ytarget);
    sourceT1CrowdDis = CrowdingDistance(sourcePopObjT1,T1ysource);
    sourceT2CrowdDis = CrowdingDistance(sourcePopObjT2,T2ysource);
    
    targetT1CrowdDis(find(isinf(targetT1CrowdDis)))=-inf;
    targetT2CrowdDis(find(isinf(targetT2CrowdDis)))=-inf;
    sourceT1CrowdDis(find(isinf(sourceT1CrowdDis)))=-inf;
    sourceT2CrowdDis(find(isinf(sourceT2CrowdDis)))=-inf;
    
    targetT1CrowdDis(find(isinf(targetT1CrowdDis)))=max(targetT1CrowdDis);
    targetT2CrowdDis(find(isinf(targetT2CrowdDis)))=max(targetT2CrowdDis);
    sourceT1CrowdDis(find(isinf(sourceT1CrowdDis)))=max(sourceT1CrowdDis);
    sourceT2CrowdDis(find(isinf(sourceT2CrowdDis)))=max(sourceT2CrowdDis);
    
    [targetT1CrowdDis,~]=mapminmax(targetT1CrowdDis,0,0.5);
    [targetT2CrowdDis,~]=mapminmax(targetT2CrowdDis,0,0.5);
    [sourceT1CrowdDis,~]=mapminmax(sourceT1CrowdDis,0,0.5);
    [sourceT2CrowdDis,~]=mapminmax(sourceT2CrowdDis,0,0.5);
    
    T1ytarget=T1ytarget-targetT1CrowdDis;
    T2ytarget=T2ytarget-targetT2CrowdDis;
    T1ysource=T1ysource-sourceT1CrowdDis;
    T2ysource=T2ysource-sourceT2CrowdDis;
    
    T1sim=corr(T1ysource',T1ytarget','type','spearman');
    T2sim=corr(T2ysource',T2ytarget','type','spearman');
    
    precentageT1=max([0,T1sim]);
    numT1=floor(length(sourcrpopulation_T1)*precentageT1);
    precentageT2=max([0,T2sim]);
    numT2=floor(length(sourcrpopulation_T2)*precentageT2);
    if pop1-numT1>0
        for i=1:pop1-numT1
            randompopulationT1(i)=Chromosome;
            randompopulationT1(i)=initialize(randompopulationT1(i),dim);
            randompopulationT1(i).skill_factor=1;
        end
        initPop=[sourcrpopulation_T1(1:numT1),randompopulationT1];
    else
        initPop=[sourcrpopulation_T1(1:numT1)];
    end
    if pop2-numT2>0
        for i=1:pop2-numT2
            randompopulationT2(i)=Chromosome;
            randompopulationT2(i)=initialize(randompopulationT2(i),dim);
            randompopulationT2(i).skill_factor=2;
        end
        initPop=[initPop,sourcrpopulation_T2(1:numT2),randompopulationT2];
    else
        initPop=[initPop,sourcrpopulation_T2(1:numT2)];
        
    end
    
    
    for i=1:pop
        population(i)=evaluate(t,initPop(i),L1,U1,mainTask,L2,U2,auxiliaryTask,dim1,dim2);
        population(i).front=[];
        population(i).rank=[];
        population(i).CD=[];
        population(i).dominationcount=0;
        population(i).dominatedset=[];
        population(i).dominatedsetlength=0;
    end
    
end

%%%%%%%%%%

population_T1=population([population.skill_factor]==1);
population_T2=population([population.skill_factor]==2);
no_of_objs_T1 = length(population_T1(1).objs_T1);
no_of_objs_T2 = length(population_T2(1).objs_T2);
[population_T1,frontnumbers]=SolutionComparison.nondominatedsort(population_T1,pop1,no_of_objs_T1);
[population_T1,~]=SolutionComparison.diversity(population_T1,frontnumbers,pop1,no_of_objs_T1);
[population_T2,frontnumbers]=SolutionComparison.nondominatedsort(population_T2,pop2,no_of_objs_T2);
[population_T2,~]=SolutionComparison.diversity(population_T2,frontnumbers,pop2,no_of_objs_T2);
population(1:pop1) = population_T1;
population(pop1+1:pop) = population_T2;

for generation=1:gen
    rndlist=randperm(pop);
    population=population(rndlist);
    for i = 1:pop 
        parent(i)=Chromosome();
        p1=1+round(rand(1)*(pop-1));
        p2=1+round(rand(1)*(pop-1));
        if population(p1).rank < population(p2).rank
            parent(i).rnvec=population(p1).rnvec;
            parent(i).skill_factor=population(p1).skill_factor;
        elseif population(p1).rank == population(p2).rank
            if rand(1) <= 0.5
                parent(i).rnvec=population(p1).rnvec;
                parent(i).skill_factor=population(p1).skill_factor;
            else
                parent(i).rnvec=population(p2).rnvec;
                parent(i).skill_factor=population(p2).skill_factor;
            end
        else
            parent(i).rnvec=population(p2).rnvec;
            parent(i).skill_factor=population(p2).skill_factor;
        end
    end
    count=1;
    for i=1:2:pop-1 % Create offspring population via mutation and crossover
        child(count)=Chromosome;
        child(count+1)=Chromosome;
        p1=i;
        p2=i+1;
        if parent(p1).skill_factor==parent(p2).skill_factor
            [child(count).rnvec,child(count+1).rnvec]=Evolve.crossover(parent(p1).rnvec,parent(p2).rnvec,muc,dim,prob_vswap);
            child(count).rnvec = Evolve.mutate(child(count).rnvec,mum,dim,1/dim);
            child(count+1).rnvec=Evolve.mutate(child(count+1).rnvec,mum,dim,1/dim);
            child(count).skill_factor=parent(p1).skill_factor;
            child(count+1).skill_factor=parent(p2).skill_factor;
        else
            if rand(1)<rmp
                [child(count).rnvec,child(count+1).rnvec]= Evolve.crossover(parent(p1).rnvec,parent(p2).rnvec,muc,dim,prob_vswap);
                child(count).rnvec = Evolve.mutate(child(count).rnvec,mum,dim,1/dim);
                child(count+1).rnvec=Evolve.mutate(child(count+1).rnvec,mum,dim,1/dim);
                child(count).skill_factor=round(rand(1))+1;
                child(count+1).skill_factor=round(rand(1))+1;
            else
                child(count).rnvec = Evolve.mutate(parent(p1).rnvec,mum,dim,1);
                child(count+1).rnvec=Evolve.mutate(parent(p2).rnvec,mum,dim,1);
                child(count).skill_factor=parent(p1).skill_factor;
                child(count+1).skill_factor=parent(p2).skill_factor;
            end
        end
        count=count+2;
    end
    for i=1:pop
        child(i)=evaluate(t,child(i),L1,U1,mainTask,L2,U2,auxiliaryTask,dim1,dim2);
    end
    population=reset(population,pop);
    intpopulation(1:pop)=population;
    intpopulation(pop+1:2*pop)=child;
    intpopulation_T1=intpopulation([intpopulation.skill_factor]==1);
    intpopulation_T2=intpopulation([intpopulation.skill_factor]==2);
    T1_pop=length(intpopulation_T1);
    T2_pop=length(intpopulation_T2);
    
    [intpopulation_T1,frontnumbers]=SolutionComparison.nondominatedsort(intpopulation_T1,T1_pop,no_of_objs_T1);
    [intpopulation_T1,~]=SolutionComparison.diversity(intpopulation_T1,frontnumbers,T1_pop,no_of_objs_T1);

    intpopulation_T1=NDS_CV(intpopulation_T1,mainTask);

    [intpopulation_T2,frontnumbers]=SolutionComparison.nondominatedsort(intpopulation_T2,T2_pop,no_of_objs_T2);
    [intpopulation_T2,~]=SolutionComparison.diversity(intpopulation_T2,frontnumbers,T2_pop,no_of_objs_T2);
    population(1:pop1) = intpopulation_T1(1:pop1);
    population(pop1+1:pop) = intpopulation_T2(1:pop2);
    
end

T1_f1=[];
T1_f2=[];
T1_f3=[];
T1_f4=[];
T2_f1=[];
T2_f2=[];
T2_f3=[];
T2_f4=[];
for i=1:pop
    if population(i).skill_factor==1
        if no_of_objs_T1 == 6
            T1_f1=[T1_f1,population(i).objs_T1(1)];
            T1_f2=[T1_f2,population(i).objs_T1(2)];
            T1_f3=[T1_f3,population(i).objs_T1(3)];
            T1_f4=[T1_f4,population(i).objs_T1(4)];
        else
            T1_f1=[T1_f1,population(i).objs_T1(1)];
            T1_f2=[T1_f2,population(i).objs_T1(2)];
        end
    else
        if no_of_objs_T2 == 4
            T2_f1=[T2_f1,population(i).objs_T2(1)];
            T2_f2=[T2_f2,population(i).objs_T2(2)];
            T2_f3=[T2_f3,population(i).objs_T2(3)];
            T2_f4=[T2_f4,population(i).objs_T2(4)];
        else
            T2_f1=[T2_f1,population(i).objs_T2(1)];
            T2_f2=[T2_f2,population(i).objs_T2(2)];
        end
        

    end
end


if mainTask.D==6
    POF=[T1_f1;T1_f2;T1_f3;T1_f4;];
    Task2POF=[T2_f1;T2_f2;T2_f3;T2_f4;];
else
    Task2POF=[T2_f1;T2_f2];
    POF=[T1_f1;T1_f2];
end

POF=[POF,Task2POF];
POF = rm_dominated(POF')';



end


function pop=NDS_CV(population_init,Problem)
temppop=[];
for i=1:length(population_init)
    temppop(i,:)=[population_init(i).rnvec,population_init(i).objs_T1];
end
V=Problem.D;
M=Problem.M;
[population front]=NDS_CD_cons(M,V,temppop);
for i=1:size(population,1)
    index=0;
    for j=1:length(population_init)
        if ismember(population(i,1:V),population_init(j).rnvec,'rows')
            index=j;
            break
        end
    end
    pop(i)=Chromosome();
    population_init=population_init';
    pop(i).rnvec=population_init(index).rnvec;
    pop(i).front=population(i,end-1);
    pop(i).objs_T1=population_init(index).objs_T1;
    pop(i).skill_factor=population_init(index).skill_factor;
    pop(i).rank=population_init(index).rank;
    pop(i).dominationcount=population_init(index).dominationcount;
    pop(i).dominatedset=population_init(index).dominatedset;
    pop(i).dominatedsetlength=population_init(index).dominatedsetlength;
end



end