clear all
load M 
load data 
for run = 1:1
    f1='ZDT4-RC';

    L1=-5*ones(1,10);
    U1=5*ones(1,10);
    L1(1)=0;U1(1)=1;
    pop1=50;
    f2='ZDT4-A';

    L2=-32*ones(1,10);
    U2=32*ones(1,10);
    L2(1)=0;U2(1)=1;
    pop2=50; 
    rmp=1; 
    gen = 250; 
    muc = 10; 
    mum = 10; 
    prob_vswap = 0; 
    store_Instance2(run).data = MOMFEA_v1(L1,U1,f1,pop1,L2,U2,f2,pop2,rmp,M,gen,muc,mum,prob_vswap,data);    
end
save('store_Instance2','store_Instance2');