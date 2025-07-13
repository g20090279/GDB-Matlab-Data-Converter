%%% Test data to verify modified code

numTests = 19;
testbench = cell(numTests,1);
newResult = cell(numTests,1);

benchData1  = load('testbench/benchData01.mat');
benchData2  = load('testbench/benchData02.mat');
benchData3  = load('testbench/benchData03.mat');
benchData4  = load('testbench/benchData04.mat');
benchData5  = load('testbench/benchData05.mat');
benchData6  = load('testbench/benchData06.mat');
benchData7  = load('testbench/benchData07.mat');
benchData8  = load('testbench/benchData08.mat');
benchData9  = load('testbench/benchData09.mat');
benchData10 = load('testbench/benchData10.mat');
benchData11 = load('testbench/benchData11.mat');
benchData12 = load('testbench/benchData12.mat');
benchData13 = load('testbench/benchData13.mat');
benchData14 = load('testbench/benchData14.mat');
benchData15 = load('testbench/benchData15.mat');
benchData16 = load('testbench/benchData16.mat');
benchData17 = load('testbench/benchData17.mat');
benchData18 = load('testbench/benchData18.mat');
benchData19 = load('testbench/benchData19.mat');

for i = 1:numTests
    eval(['testbench{',num2str(i),'}  = benchData',num2str(i),'.logData']);
end

%%% Convert File
convGdbLog2Mat('testbench/gdb.log.testbench.txt',1);

%%% Compare Results
if exist('gdb_log_testbench.m', 'file') == 2
    run('gdb_log_testbench.m');
    convData1 = var1;
    convData2 = var2;
    convData3 = var3;
    convData4 = var4;
    convData5 = var5;
    convData6 = var6;
    convData7 = var7;
    convData8 = var8;
    convData9 = var9;
    convData10 = var10;
    convData11 = var11;
    convData12 = var12;
    convData13 = var13;
    convData14 = var14;
    convData15 = var15;
    convData16 = var16;
    convData17 = var17;
    convData18 = var18;
    convData19 = var19;
else
    % logNewAgcGain = load('.mat');
    convData1 = load('gdb.log.testbench_var1.mat').logData;
    convData2 = load('gdb.log.testbench_var2.mat').logData;
    convData3 = load('gdb.log.testbench_var3.mat').logData;
    convData4 = load('gdb.log.testbench_var4.mat').logData;
    convData5 = load('gdb.log.testbench_var5.mat').logData;
    convData6 = load('gdb.log.testbench_var6.mat').logData;
    convData7 = load('gdb.log.testbench_var7.mat').logData;
    convData8 = load('gdb.log.testbench_var8.mat').logData;
    convData9= load('gdb.log.testbench_var9.mat').logData;
    convData10 = load('gdb.log.testbench_var10.mat').logData;
    convData11 = load('gdb.log.testbench_var11.mat').logData;
    convData12 = load('gdb.log.testbench_var12.mat').logData;
    convData13 = load('gdb.log.testbench_var13.mat').logData;
    convData14 = load('gdb.log.testbench_var14.mat').logData;
    convData15 = load('gdb.log.testbench_var15.mat').logData;
    convData16 = load('gdb.log.testbench_var16.mat').logData;
    convData17 = load('gdb.log.testbench_var17.mat').logData;
    convData18 = load('gdb.log.testbench_var18.mat').logData;
    convData19 = load('gdb.log.testbench_var19.mat').logData;
end

for i = 1:numTests
    if eval(['isequal(testbench{', num2str(i), '}, convData', num2str(i), ')'])
        disp(['Test ', num2str(i), ' passed.']);
    else
        disp(['Test ', num2str(i), ' failed.']);
    end
end
