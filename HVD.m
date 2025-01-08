function Score = HVD(PopObj,PF)
%% PF is the true POF
    r = max(PF)+0.5;
    Score = hypervolume(PF,r) - hypervolume(PopObj,r);
end