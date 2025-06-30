function [data, time] = load_data()
    csf_data_file1 = 'data/blattner_wake_conc.csv';
    csf_data_file2 = 'data/lucey_wake_conc.csv';
    csf_data_file3 = 'data/liu_csf_wake_conc.csv';
    plasma_data_file1 = 'data/liu_plasma_wake_conc.csv';

    csf_data1 = readtable(csf_data_file1);
    csf_data2 = readtable(csf_data_file2);
    csf_data3 = readtable(csf_data_file3);
    plasma_data1 = readtable(plasma_data_file1);

    time.exp1 = csf_data1.Time;
    time.exp2 = csf_data2.Time;
    time.exp3 = csf_data3.Time;

    % Extracting conc
    data.csf_conc_exp1 = csf_data1.Concentration;
    data.csf_conc_exp2 = csf_data2.Concentration;
    data.csf_conc_exp3 = csf_data3.Concentration;
    data.plasma_conc_exp1 = plasma_data1.Concentration;

    % Extracting error
    data.csf_lsd1 = csf_data1.LSD;
    data.csf_lsd2 = csf_data2.LSD;
    data.csf_lsd3 = csf_data3.LSD;
    data.csf_usd1 = csf_data1.USD;
    data.csf_usd2 = csf_data2.USD;
    data.csf_usd3 = csf_data3.USD;
    data.plasma_lsd1 = plasma_data1.LSD;
    data.plasma_usd1 = plasma_data1.USD;

    % Calculating SD
    data.csf_std1 = (data.csf_usd1 - data.csf_lsd1)/2;
    data.csf_std2 = (data.csf_usd2 - data.csf_lsd2)/2;
    data.csf_std3 = (data.csf_usd3 - data.csf_lsd3)/2;
    data.plasma_std1 = (data.plasma_usd1 - data.plasma_lsd1)/2;
end
