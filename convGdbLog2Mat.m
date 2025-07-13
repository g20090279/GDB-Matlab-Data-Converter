function convGdbLog2Mat (fileName, outFileMode, inputMode)
% CONVGDBLOG2MAT converts logging data from GDB into .mat or .m file. It 
% works with C++ Eigen library, where the data looks like
%          Vector < Vector < Vector < ... Eigen::Matrix<> ... > > >,
% i.e. only one Eigen::Matrix is capsulated by multiple std::vector. Note
% that the dimension of Eigen::Matrix in each std::vector can be different.
%
% This function searches only for the files staring with 'gdb.log.' if the
% input argument 'fileName' is not specified.
%
% Inputs:
%   - fileName: (optional) A char array specifying the file.
%   - outFileMode: the format of output data
%       * 0: one .mat file for each data
%       * 1: one .m file for all interface data
%   - inputMode: which data types are processed
%       * 0: convert only Eigen::Matrix data
%       * 1: convert only all types of data
%
% Outputs: (in the .mat file)
%   - logData: The converted data with same dimension
%   - dimStdVec: A column vector indicates number of elements for all
%       std::vector.
%
%       Example:
%
%       dimStdVec = [14; 4; 1] for
%       Vector<Vector<Vector<Eigen::Matrix<> ...>>>
%         ^      ^      ^
%        14      4      1     <- the length of std::vector
%   - dimEigMatInVec: A two-column matrix, where each row indicates the
%       dimension of Eigen::Matrix of the corresponding std::vector element.
%       Therefore, The number of rows equals to the product of all elements
%       in dimStdVec. Note that the row index starts from the most inner
%       std::vector.
%
%       Example:
%
%       dimStdVec = [1,2,3] for
%       Vector<Vector<Vector<Eigen::Matrix<> ...>>>
%         ^      ^      ^
%         1      2      3     <- the length of std::vector
%
%       dimEigMatInVec = [
%           12, 2     --
%           12, 2      |--> dimStdVec(1,1,:)
%           12, 2     --
%           14, 1     --
%           14, 1      |--> dimStdVec(1,2,:)
%           14, 1     --
%       ]
%
% Examples:
%   - convGdbLog2Mat('~/data/mylog.txt') will convert the file '~/data/mylog.txt'
%   - convGdbLog2Mat('~/data/') will convert all files under the directory '~/data/' with name prefix 'gdb.log.'.
%   - Put the script in the same folder as log files, and run the script without any input. This will also convert all the files in the same directory.

if nargin < 3, inputMode = 0; end
if nargin < 2, outFileMode = 0; end

%%% temporary check
if inputMode>0, error('inputMode not supported!'); end
if outFileMode>1, error('outFileMode not supported'); end

logFile = {};   % input log files, can be multiple
out = '';

if nargin > 0 && ~isempty(fileName) % check if a specific file given, it can has any file name
    if isfile(fileName)
        logFile = {fileName};
    end
else
    fileName = '.';
end

if exist(fileName,'dir') == 7  % 7 = directory, search all files with prefix gdb.log.
    fileSearch = [fileName,'/gdb.log.*'];
    disp(['Searching for ',fileSearch]);
    fileInfo = dir(fileSearch);
    cFiles = 0;
    if ~isempty(fileInfo)
        for iFile=1:length(fileInfo)
            % Only Process File with Name Beginning with gdb.log and Not with
            % File Extension .mat
            [~,fileName,fileExt] = fileparts(fileInfo(iFile).name);
            if ~strcmp(fileExt,'.mat') && strcmp(fileName(1:8),'gdb.log.')
                cFiles = cFiles + 1; % counter
                fileFullPath = [fileInfo(iFile).folder,'\' fileInfo(iFile).name];
                if isunix
                    fileFullPath(strfind(fileFullPath,'\')) = '/';
                end
                logFile = [logFile, fileFullPath];
                continue;
            end
        end
    end
end

disp([num2str(length(logFile)),' File(s) found.'])

if ~isempty(logFile)

    % Begin Converting
    for iFile = 1:length(logFile)

        disp(['Processing File ', logFile{iFile}, ' ...']);
        [~,convertingFileName,~] = fileparts(logFile{iFile});
        fid = fopen(logFile{iFile});
        tLine = fgetl(fid);

        % Initialize values
        isFound = false;      % if a gdb log file for dump data from GDB is found
        isOneline = false;    % if Eigen data is one-line data or spans into multiple lines
        isComplete = false;   % if the data is complete (need for multiple line)
        isFirstWrite = true; % use for outFileMode==1, first write to empty file

        % Maintain a 3-line buffer
        % Example log setting file:                                       Example log file:
        % | ...                                                           | ...                             |
        % | ...                                                           | +print "IF:in_var"              |
        % | print "IF:in_var         "         ------------------>        | $4 = "IF:in_var"                | <- tLinePrev2
        % | print var                               (log file)            | +print var                      | <- tLinePrev1
        % | ...                                                           | $5 = Eigen::Matrix<...          | <- tLine
        % | ...                                                           | ...                             |
        %
        % Print the IF name before the variable itself with the indicator "IF:"
        tLinePrev2 = '';
        tLinePrev1 = '';

        while ischar(tLine)
            % Note 1: the data structure is always: the innermost is Eigen, capsulated by multiple std::vector.
            % Note 2: the eigen matrix has maximum two dimensions.
            % Note 3: matrix element listing begins from innermost.

            %%% Detect data
            % - Eigen data: .m & .mat
            % - scalar number (integer, double, complex, positive/negative): .m
            % - string (boolean, enumeration): .m
            if ~isFound  % searching for first line

                %%% Detect pattern example: $2 = std::vector of length 14, capacity 14 = { std::vector of length 5, capacity 4 = {Eigen::Matrix<std::complex<double>,1200,2,ColMajor> (data ptr: 0x4fcfe80) = {[0,0] = {_M_value = -0.063488650798086924 + 0.10254959057443244 * I}, [1,0] = ...
                %                           \-/   \--------------------------------------v------------------------------------/   \-------------------------v------------------------/
                %                       1st token                                    2nd token                                                          3rd & 4th token
                patternEigen = regexp(tLine, '^[\w~"]*\$(\d+)\s*=\s*(std::vector of length\s*\d+[,\s]*capacity\s+\d+\s*=\s*{)*Eigen::Matrix<[\w<>:]+,(\d+),(\d+),\w+>', 'tokens', 'once');

                %%% If (1) no pattern is found, or (2) data is not starting with "[~"]$numbers =', or (3) no Eigen data (identified by valid dimensions)
                if isempty(patternEigen) || isempty(patternEigen{1}) || isempty(patternEigen{3}) || isempty(patternEigen{4})
                    processNextLine;
                    continue;
                end

                %%% Check if the data is one line or spans multiple lines
                numLeftBracket = sum(tLine=='{');
                numRightBracket = sum(tLine=='}');
                if numLeftBracket == numRightBracket
                    isOneline = true;
                end

                %%% Detect if there are std::vector(s) encapsulating Eigen::Matrix
                resStdVec = regexp(patternEigen{2}, 'std::vector of length ([\d]+), capacity ([\d]+)', 'tokens');
                numStdVec = length(resStdVec);  % Number of std::vectors wrapping Eigen::Matrix
                dimStdVec = zeros(numStdVec,1);
                for i = 1:numStdVec
                    dimStdVec(i) = str2double(resStdVec{i}{1});
                end
                % Note: if dimStdVec is empty, numElStdVec is 1, which is correct in our scenario.
                numElStdVec = prod(dimStdVec);

                % Record the dimensions of Eigen::Matrix in each std::vector, which may be different
                dimEigMatInVec = zeros(numElStdVec, 2);

                if isOneline % Has enough info by using this line
                    eigMatAll = regexp(tLine, 'Eigen::Matrix<[\w:<>]+,(\d+),(\d+),\w+>','tokens');  % Record all Eigen::Matrix

                    if length(eigMatAll) ~= numElStdVec  % #occurrence of Eigen::Matrix should equal to #elements of all std::vector
                        disp('Error: Occurrence of Eigen::Matrix is not equal the sum of std::vector elements!');
                        tLine = fgetl(fid);
                        continue;
                    end

                    for i = 1:numElStdVec  % Put from cell to vector
                        dimEigMatInVec(i,1) = str2double(eigMatAll{i}{1});
                        dimEigMatInVec(i,2) = str2double(eigMatAll{i}{2});
                    end
                else  % Need to go next lines until this data block ends or new data block begins or end of file
                    % Save the current position
                    posFid = ftell(fid);

                    % Use the number of '{' and '}' to check if the current data block ends
                    cnt = 0;   % counter for multi-line case
                    isNewData = false;
                    tLineTmp = tLine;

                    while  numLeftBracket~=numRightBracket && ~isNewData && cnt<numElStdVec
                        % Search only for Eigen::Matrix
                        eigMatAll = regexp(tLine, 'Eigen::Matrix<[\w:<>]+,(\d+),(\d+),[\w]+>','tokens');  % Record all dim of Eigen::Matrix

                        % Get new line and check "{" and "}"
                        tLine = fgetl(fid);
                        numLeftBracket = numLeftBracket + sum(tLine=='{');
                        numRightBracket = numRightBracket + sum(tLine=='}');

                        if isempty(eigMatAll)
                            continue;
                        end

                        if length(eigMatAll) ~= 1
                            disp('Multi-line data can only have one Eigen::Matrix at one line.');
                            break;
                        end

                        % Record all dim of Eigen::Matrix
                        cnt = cnt + 1;
                        dimEigMatInVec(cnt,1) = str2double(eigMatAll{1}{1});
                        dimEigMatInVec(cnt,2) = str2double(eigMatAll{1}{2});

                        % Check data block ending indicators. If the current data block ends, break the while.
                        if ~ischar(tLine), isNewData = true; continue; end

                        findTextGdb = regexp(tLine, 'gdb = std::vector of length','tokens', 'once');
                        if ~isempty(findTextGdb), isNewData = true; continue; end

                        findTextDollar = regexp(tLine, '\$(\d+) = std::vector of length','tokens', 'once');
                        if ~isempty(findTextDollar)
                            if str2double(findTextDollar{1}) ~= indCommand
                                isNewData = true;
                            end
                        end
                    end

                    if cnt ~= numElStdVec
                        disp('Error: Occurrence of Eigen::Matrix is not equal the sum of std::vector elements!');
                        continue;
                    else  % Data is fully collected
                        % Go back to first line of the data
                        fseek(fid,posFid,'bof');
                        tLine = tLineTmp;
                    end
                end

                dimEigMatMax = max(dimEigMatInVec,[],1);  % maximum dimension of eigen matrix

                % Dimension of the data. Eigen::Matrix column-weise. Eigen data is the innermost
                dimAll = [dimStdVec.', fliplr(dimEigMatMax)];   % use this format vector(outer)..vector(inner),Eigen::Matrix.numCols,Eigen::Matrix.numRows is to preserve column dimension by preventing auto-squeezing in Matlab
                numDimAll = length(dimAll);

                % Initialize Variables
                logData     = zeros(prod(dimAll),1);
                isFound     = true;
                isComplete  = false;
                isError     = false;
                errMsg      = '';

                % Use following indices to check if enough data are collected for each column of a Eigen::Matrix.
                % Align the writing index iData when the dimension of a Eigen::Matrix is smaller than the maximum dimension (dimEigMat)
                cntData     = 1;      % index writing to the result matrix
                cntVec      = 1;      % index of std::vector, the dimension of Eigen::Matrix in this std::vector may different
                currColInd  = 1;      % current column index inside one Eigen::Matrix
                currRowInd  = 1;      % current row index inside one Eigen::Matrix
                prevRowInd  = 1;      % previous row index inside one Eigen::Matrix
                iData       = 1;      % index of all Eigen::Matrix

                if isFound
                    % Get variable name information from GDB log
                    indCommand = str2double(patternEigen{1});
                    varNameIf = regexp(tLinePrev2, '^[~"]*\$\d+\s*=[\s]*[\\"]*IF\s*:\s*(\w+)[\\"].*', 'tokens', 'once');
                    varNameLog = regexp(tLinePrev1, 'p(?:rint)?\s*(?:[\w]+[->.:]+)*([\w]+)', 'tokens', 'once');
                    if ~isempty(varNameIf)    % When IF name is given, this will overwrite the variable name
                        convertingDataName = ['var', num2str(indCommand), '_', varNameIf{1}];
                    else        
                        if ~isempty(varNameLog)
                            convertingDataName = ['var', num2str(indCommand), '_', varNameLog{1}];
                        else
                            convertingDataName = ['var', num2str(indCommand)];
                        end
                    end
                    disp(['Data ', convertingDataName, ' is found! Start converting ...']);
                end
            end

            if isFound  % If data block is found

                if isOneline % No need to read next lines
                    %%% Extract data
                    % normal integer format: [0] = -5
                    % normal double format: [0,0] = -19.952623149688797
                    % normal complex-valued double format: [0,0] = {_M_value = 0.14589548774745745 + -0.13207671786275996i}
                    % fxp format: [0] = {_M_real = -1728, _M_imag = 1312}
                    resValue = regexp(tLine, '\[(\d+)(?:,)?(\d+)*\]\s*=\s*(?:{(_M_value|_M_real)\s*=\s*)?(-?\s*\d*\.?\d*e?-?\d*)(?:\s*(\+|,?\s*_M_imag\s*=*)\s*)?(-?\s*\d*\.?\d*e?-?\d*)*', 'tokens');

                    if length(resValue) ~= sum(prod(dimEigMatInVec,2))   % If didn't capture all data
                        isFound = false;
                        disp('Error: not enough data!');
                    end

                    for iData = 1:length(resValue)  % Write data
                        % If this column length is not equal to max column length, jump to max column length
                        checkDataIndex;
                        if isError
                            disp(errMsg);
                            break;
                        end

                        % Write Data, complex-valued or real-valued
                        writeData;
                        cntData = cntData + 1;

                        % Check if enough data for this Eigen::Matrix column, and if is complete
                        checkEnoughDataThisColumn;
                    end
                else
                    % Check if out of data block
                    findTextDollar = regexp(tLine, '\$([\d]+) = std::vector of length','tokens', 'once');
                    if ~isempty(findTextDollar)
                        if str2double(findTextDollar{1}) ~= indCommand
                            isFound = false;
                            disp('Error: data block terminates.');
                        end
                    end

                    findTextGdb = regexp(tLine, 'gdb = std::vector of length', 'once');
                    if ~isempty(findTextGdb)
                        isFound = false;
                        disp('Error: data block terminates.');
                    end

                    % Capture Eigen::Matrix data.
                    % Pattern 1 - vector with real valued data:    Eigen::Matrix<int,864,1,ColMajor> (data ptr: 0x25276f0) = {[0] = 3"                                                                                % for one-column data, [0] instead of [0,0]
                    % Pattern 2 - vector with complex-valued data: Eigen::Matrix<std::complex<double>,144,1,ColMajor> (data ptr: 0x1d5a2f0) = {[0] = {_M_value = 0.44140625 + -0.8759765625 * I"                      % _M_value is present for complex-valued data
                    % Pattern 3 - matrix with complex-valued data: Eigen::Matrix<std::complex<double>,288,2,ColMajor> (data ptr: 0x1d927c0) = {[0,0] = {_M_value = 0.21752604580125023 + -0.43954164224987846 * I"
                    % Note 1: For column data, the second element of resValue is empty.
                    % Note 2: For real-valued data, the fourth element of resValue is empty.
                    resValue = regexp(tLine, '\[(\d+)(?:,)?(\d+)*\]\s*=\s*(?:{_M_value\s*=\s*)?(-?\s*\d*\.?\d*e?-?\d*)(?:\s*\+\s*)?(-?\s*\d*\.?\d*e?-?\d*)*', 'tokens');

                    if ~isempty(resValue)
                        if length(resValue) ~= 1  % For multi-line data, only one element is present on each line.
                            isFound = false;
                            disp('Error: For multi-line data, only one data appears on each line!');
                            continue;
                        end

                        checkDataIndex;
                        if isError
                            disp(errMsg);
                            continue;
                        end

                        % Write data, complex-valued or real-valued
                        writeData
                        cntData = cntData + 1;

                        % Check if enough data for this Eigen::Matrix column, and if is complete
                        checkEnoughDataThisColumn
                    end
                end

                %%% If data is complete in this data block
                if isComplete
                    % Reshape data back to origin dimensions
                    if numDimAll > 1
                        % Note: reshape will sequeeze out one-dim
                        logData = reshape(logData, fliplr(dimAll));
                        logData = permute(logData, [length(size(logData)):-1:3,1,2]);
                        logData = reshape(logData, [dimStdVec.',dimEigMatMax]);
                    end

                    %%% Output this data, 3 data for each Eigen data
                    % 1. logData - the data
                    % 2. dimStdVec - the number of std::vector encapulating Eigen
                    % 3. dimEigMatInVec - the length of Eigen in each std::vector, this can vary in each std::vector

                    disp('Conversion done successfully!');

                    if outFileMode == 0
                        if ~isempty(varNameIf)  % if IF name is given, no file name prefix. Warning: data will be overwritten.
                            save([convertingDataName,'.mat'],'logData','dimStdVec','dimEigMatInVec');
                        else
                            save([convertingFileName, '_', convertingDataName,'.mat'],'logData','dimStdVec','dimEigMatInVec');  % if IF name is not given, 
                        end
                    elseif outFileMode == 1
                        % Replace with '_' to have valid file name
                        tmpFileName = convertingFileName;
                        tmpFileName = [regexprep(tmpFileName, '[.\s]+', '_'), '.m'];

                        % If first write, make the file empty. Otherwise, append.
                        if isFirstWrite
                            fidw = fopen(tmpFileName,'w');
                            isFirstWrite = false;
                            fclose(fidw);
                        end
                        writeEigen2DotmFile(tmpFileName, convertingDataName, logData, dimStdVec, dimEigMatInVec);
                    else
                        disp('The output file format option is not supported');
                    end

                    % Reset parameters
                    isFound     = false;
                    isComplete  = false;
                    isOneline   = false;
                end
            end

            % Read the new line
            tLine = fgetl(fid);

        end

        % Finish converting this file
        fclose(fid);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%
%   Nested Functions   %
%%%%%%%%%%%%%%%%%%%%%%%%

function processNextLine
    % Save current two lines, read the new line
    tLinePrev2 = tLinePrev1;
    tLinePrev1 = tLine;
    tLine = fgetl(fid);
end

function checkDataIndex
    if isempty(resValue{iData}{2})
        currColInd = 1;
    else
        currColInd = str2double(resValue{iData}{2}) + 1;
    end

    prevRowInd = currRowInd;
    currRowInd = str2double(resValue{iData}{1}) + 1;

    if cntData ~= (cntVec-1)*prod(dimEigMatMax) + (currColInd-1)*dimEigMatMax(1) + currRowInd % Track global counter and the corresponding index in a column
        isFound = false;
        isError = true;
        errMsg = 'Error: wrong data index!';
    end

    if prevRowInd > currRowInd && prevRowInd ~= dimEigMatMax(1)    % New column, check if lack of Eigen::Matrix data of last column
        isFound = false;
        isError = true;
        errMsg = ['Error: The ', num2str(cntVec), ' std::vector has not enough data!'];
    end
end

function writeData
    if isempty(resValue{iData}{4}) % Real-valued
        logData(cntData) = str2double(resValue{iData}{3});
    else % Complex-valued
        logData(cntData) = str2double(resValue{iData}{3}) + 1i * str2double(resValue{iData}{4});
    end
end

function checkEnoughDataThisColumn
    if mod(cntData-1, dimEigMatMax(1)) == dimEigMatInVec(cntVec,1) || mod(cntData-1, dimEigMatMax(1)) == 0   % Enough data for this column
        cntData     = cntData + dimEigMatMax(1) - dimEigMatInVec(cntVec,1);   % jump to the max row dim
        currRowInd  = dimEigMatMax(1);
        if currColInd == dimEigMatInVec(cntVec,2)  % If the last column, jump to the max col dim
            cntData = cntData + dimEigMatMax(1) * ( dimEigMatMax(2) - dimEigMatInVec(cntVec,2) );
            cntVec  = cntVec + 1;   % std:vector index increases 1
        end
        % Check if all data are captured
        if cntVec == numElStdVec + 1
            isComplete = true;
        end
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Local Subfunctions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeEigen2DotmFile(filename, varname, data, dimVec, dimEigInVec)
    fid = fopen(filename, 'a');

    % Print the length of std::vector from outer to inner, as comment code
    fprintf(fid, '\n%% For data %s, number of std:vector is %d', varname, length(dimVec));
    if isempty(dimVec)
        fprintf(fid, '.');
    else
        fprintf(fid, ', the length of each std:vector from outer to inner is [');
        for i = 1:length(dimVec)-1
            fprintf(fid, '%d,', dimVec(i));
        end
        fprintf(fid, '%d].', dimVec(end));
    end

    % Print the dimension of Eigen Matrix in each std::vector element, from inner to outer
    fprintf(fid, '\n%% For data %s, the dimension of Eigen Matrix for each std::vector element (from inner to outer) is: [', varname);
    for i = 1:size(dimEigInVec,1)-1
        fprintf(fid, '%d,%d;', dimEigInVec(i,1), dimEigInVec(i,2));
    end
    fprintf(fid, '%d,%d]\n', dimEigInVec(end,1), dimEigInVec(end,2));

    % Print data
    % fprintf(fid, '%s(i,:)=[%.17f,', varname, data(i,j,:));
    dim = size(data);
    numDim = length(dim);

    if numDim > 2  % multiple rows in .m file
        numItem = prod(dim(1:end-2));
        dimExtra = zeros(numItem, numDim-2);  % every row write the data from last dimension

        % Loop for each dimension larger than 2
        for i = 1:numDim-2
             tmpDim = repmat(1:dim(i), prod(dim(i+1:end-2)), 1);
             dimExtra(:,i) = repmat(tmpDim(:), prod(dim(1:i-1)), 1);
        end

        % Print the in matrix form
        for i = 1:numItem
            % The first row, i.e. 1st element of the second last dimension
            tmpIdx = regexprep(num2str(dimExtra(i,:)),'\s+',',');
            eval(['fprintf(fid, ''', varname, '(', tmpIdx, ',:,:)=['');']);

            for j = 1:dim(end-1) % for each row in the matrix
                % Need to handle real- or complex-valued data
                if isreal(data)
                    eval(['fprintf(fid, ''%.17g, '', data(', tmpIdx, ',', num2str(j), ',1:end));']);
                else
                    eval(['tmpData=squeeze(data(', tmpIdx, ',', num2str(j), ',1:end));']);
                    tmpData = [real(tmpData).'; imag(tmpData).'];
                    fprintf(fid, '%.17g + %.17gi, ', tmpData(:));
                end
                fprintf(fid, '; ');
            end
            fprintf(fid, '];\n');
        end

    else  % only one row in .m file
        fprintf(fid, '%s=[', varname);

        for i = 1:dim(1)
            if isreal(data)
                fprintf(fid, '%.17g, ', data(i,:));
            else
                tmpData = [real(data(i,:));imag(data(i,:))];
                fprintf(fid, '%.17g + %.17gi, ', tmpData(:));
            end
            fprintf(fid, '; ');
        end
        fprintf(fid, '];\n');
    end

    % Close file pointer
    fclose(fid);
end
