if(~exist('LeadfieldMatrix_adult.mat','file'))
    load('LeadfieldMatrix_adult_part1');
    load('LeadfieldMatrix_adult_part2');
    load('LeadfieldMatrix_adult_part3');
    load('LeadfieldMatrix_adult_part4');
    Gain = [Gain1 Gain2 Gain3 Gain4];
    Brain = [Brain1; Brain2; Brain3; Brain4];
    clear Gain1 Gain2 Gain3 Gain4
    clear Brain1 Brain2 Brain3 Brain4
    save('Lead Field Matrix Data\LeadfieldMatrix_adult.mat','Brain','Gain')
end