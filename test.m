% test for runtesst.m

findn = 10;
badseeds = [];
for seed = 2000:3000
    runtest;
    difference = eigenlist(end) - disp_list(end)
    seed
    if difference > 0.1
        badseeds = [badseeds, seed];
    end
end

