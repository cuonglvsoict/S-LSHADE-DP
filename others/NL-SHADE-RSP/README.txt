#Authors: Vladimir Stanovov, Shakhnaz  Akhmedova, Eugene Semenkin

A short readme on how to compile and run the code of NL-SHADE-RSP for CEC 2021 benchmark functions.
CEC 2021 competition website:
http://www.ntu.edu.sg/home/epnsugan/index_files/CEC2017/CEC2017.htm

All the code of the NL-SHADE-RSP is contained in "nl_shade_rsp.cpp" file.
For compilation you will also need "cec21_test_func.cpp", provided by the competition organizers.
To run the code, you will need the "input_data" folder, containing the functions shift and rotation data, as well as random seeds.

The basic compilation is simple using gcc/g++:

g++ -std=c++11 -O3 nl_shade_rsp.cpp -o nl_shade_rsp.exe

You may try to decrease execution time by adding hardware-specific commands:

g++ -std=c++11 -O3 -march=corei7-avx nl_shade_rsp.cpp -o nl_shade_rsp.exe

Please note that the compilation requires support of C++11 standard.
This will create nl_shade_rsp.exe, available for running.
The program will first run time complexity tests on your machine, writing results into "time_complexity.txt"
You may edit the .cpp file to switch it off.
Next, the main optimization loop will be started, writing data to "NL-SHADE-RSP_(###)_(F)_(DIM)D.txt", where ###, F and DIM are the benchmark code, function number and problem dimension.
As long as the the random seeds are fixed, the results should always be the same.

Should you face any problems, or have questions, my e-mail is vladimirstanovov@yandex.ru.
