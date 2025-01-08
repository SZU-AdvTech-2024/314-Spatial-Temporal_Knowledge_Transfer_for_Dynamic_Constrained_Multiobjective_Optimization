function con=configure()
con.T_parameter = [
    1 20 600
       10  20 600
      10  10  300    
                    ];%% time parameters 
con.TestFunctions = {'S1T1','S1T3','S2T1','S2T2','S2T3','S2T4','S3T1','S3T2','S3T3','S3T4','S4T1','S4T2','S4T3','S4T4'};%DCF1-14

con.popSize=100;
con.repeat=10;
con.dec=10;

end