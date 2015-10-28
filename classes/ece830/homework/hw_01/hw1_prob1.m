clear;
close all;

data_files = {'brain_data1.mat', 'brain_data2.mat'};

for i=1:2
    load(data_files{i})

    p_xA = sum(xA == xB) / length(xA);
    p_xC = sum(xC == xB) / length(xA);

    disp('P(x_A, x_C | x_B')
    disp(p_xA * p_xC)

end

