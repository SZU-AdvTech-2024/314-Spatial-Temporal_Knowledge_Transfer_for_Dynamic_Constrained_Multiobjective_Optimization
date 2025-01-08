classdef Chromosome
    
    properties
        rnvec;
        objs_T1;
        objs_T2;
        convio;
        skill_factor;
        front;
        CD;
        rank;
        dominationcount=0;
        dominatedset=[];
        dominatedsetlength=0;
    end
    
    methods
        
        function object=initialize(object,dim)
            object.rnvec=rand(1,dim);
        end
        
        function object=evaluate(t,object,L1,U1,f1,L2,U2,f2,dim1,dim2)
            M=f1.M;
            D=f2.D;
            if object.skill_factor==1
                xtemp=object.rnvec(1:dim1);                
                x=L1+xtemp.*(U1-L1);
                object.objs_T1=f1.getObj(x,M,D,t); 
                if dim1==6 
                    CV=f1.getCV(object.objs_T1,t,x);
                else    
                     CV=f1.getCV(object.objs_T1,t);
                end
               
                object.objs_T1=[object.objs_T1,CV*5];
           
            else
                xtemp=object.rnvec(1:dim2);                
                x=L2+xtemp.*(U2-L2);                 
                object.objs_T2=f2.getObj(x,M,D,t,f1);  
            end            
        end   
        
        function population=reset(population,pop)
            for i=1:pop
                population(i).dominationcount=0;
                population(i).dominatedset=[];
                population(i).dominatedsetlength=0;
            end
        end     
    end    
    
end

