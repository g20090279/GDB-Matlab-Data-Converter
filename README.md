---
Topic: GDB Logging Data to Matlab Data Converter
Description: To convert the GDB logged Eigen::Matrix Data to Matlab .mat or .m file.
Last Revised: July 13, 2025
---

# Background

Matlab mex interface enables C/C++ code in Matlab. To verify that the mex function is correctly implemented, we can log the C/C++ data from GDB debugger, and then use them in Matlab for verification. One option is to write `cout` function to print the data into files. This needs to recompile C/C++ project every time when a modification is added. In addition, the changes for printing data will be lost when the C/C++ project is updated.

To overcome these disadvantages, we can use the debugger (can be GDB, LLVM, but here refers only to GDB) to log the data into files. Moreover, the printing commands can be put into a single file, which can be used later to have full automation for data acquiring.

# Generate Log Data from C/C++ Project

1. Compile C/C++ project in debug mode.
2. Add [pretty printer plugin](https://gitlab.com/libeigen/eigen/-/tree/master/debug/gdb?ref_type=heads) (both for STL vector and Eigen) for GDB.
3. Run the program (with debug mode) until the program stops at the breakpoint. Afterwards, use the following commands to log the data into a file.

```
(gdb) set print elements 0
(gdb) set trace-commands on
(gdb) set pagination off
(gdb) set logging file gdb-log/gdb.log.mydata.txt
(gdb) set logging enabled on
(gdb) print "IF:wanted_name_of_this_variable"
(gdb) print this_variable
(gdb) print that_variable
(gdb) set logging enabled off
```

Note also that in VSCode, a prefix `-exec` before these commands has to be added.

## Automatic Generation

If you want to extract same data with different setting files, probably you don't want to enter these commands again and again. You can set up a GDB script to run GDB commands automatically. GDB can run a script with `-x` option. For example,

```bash
$ gdb -x extract_data.gdb
```

We can also use VSCode by changing the setting in `launch.json` accordingly.

# Current Limitation

1. Only the works with C++ Eigen library. The script checks the key word "Eigen::Matrix" to identify a data block. Note that the data has a type structure `Vector<Vector<Vector<...Eigen::Matrix<>...>>>`.
1. You can specify the file name in input argument. If not, it only searches for the files staring with `gdb.log`. Therefore, you should set the log file name in GDB as `gdb.log.<variablename>[.txt]`.

# Use Example

1. `convGdbLog2Mat('~/data/mylog.txt')` will convert the first data block in the file `~/data/mylog.txt` into `~/data/mylog.mat`
2. `convGdbLog2Mat('~/data/')` will convert all files (first data block) under the directory `~/data/` with name prefix `gdb.log.`.
3. Put the script in the same folder as log files, and run the script without any input. This will also convert all the files in the same directory.
4. `convGdbLog2Mat('~/data/mylog.txt',1)` will convert data to a `.m` file.

# Testbench

A testbench is ready to verify the changes in the code.

# CHANGELOG

- v0.7 **Useful feature**
  - need to clear repo history,
  - converted data can be store in a `.m` file besides `.mat` file,
  - use `IF:` flag to name the variable,
  - improve regexp to be more general to identify the data with dummy characters from GDB

- v0.6
  - fix bug for case with different Eigen::Matrix dimensions,
  - add test case for the above case.

- v0.5
  - fix error in data index checking when different size of Eigen::Matrix in std::vector element. Note that this data index checking scheme is to verify that enough data are captured for each column of Eigen::Matrix,
  - improve variable naming when prefix is present,
  - add testbench.

- v0.4 **IMPORTANT UPDATE**
  - [**IMPORTANT!**] fix error: in the data, the 'e' representing exponent is not captured,
  - fix error in capturing the dimensions of std::vector,
  - stop throwing error handles.

- v0.3
  - fix error when Eigen::Matrix is not capsulated by std::vector.

- v0.2
  - supports multi-line data and one-line data,
  - supports varying length of Eigen::Matrix.

- v0.1 initial version
  - functionality: generate `.mat` data for MATLAB from the GDB log data for Eigen library,
  - limitation: same data length for all dimensions, only first data block of each log file will be processed, work only with Eigen library.
