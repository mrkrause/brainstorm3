function varargout = process_move_raw( varargin )
% PROCESS_MOVE_RAW: Move raw files to new paths
%
% USAGE:  process_move_raw('Compute', filename, oldPath, newPath)

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2017 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Martin Cousineau, 2017

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Move raw files';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'File';
    sProcess.Index       = 1025;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    % Option: Epoched/continuous
    sProcess.options.oldpath.Comment = 'Previous parent path of raw files:';
    sProcess.options.oldpath.Type    = 'text';
    sProcess.options.oldpath.Value   = '';
    sProcess.options.newpath.Comment = 'New parent path of raw files:';
    sProcess.options.newpath.Type    = 'text';
    sProcess.options.newpath.Value   = '';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    % Options
    oldPath = sProcess.options.oldpath.Value;
    newPath = sProcess.options.newpath.Value;
    
    if isempty(oldPath) || isempty(newPath)
        bst_report('Error', sProcess, sInputs, 'Old and new parent paths are mandatory.');
        return;
    end
    if exist(newPath, 'dir') ~= 7
        bst_report('Error', sProcess, sInputs, 'The new parent path does not exist.');
        return;
    end
    
    % Convert all the files in input
    for i = 1:length(sInputs)
        % Load file
        DataFile = file_fullpath(sInputs(i).FileName);
        DataMat = in_bst_data(DataFile);
        sFile = DataMat.F;

        % Convert
        [sFile, warnMessage] = Compute(sFile, oldPath, newPath);
        
        % Error handling
        if ~isempty(warnMessage)
            bst_report('Warning', sProcess, sInputs(i), [warnMessage ' Skipping this file.']);
        else
            % Add history field
            DataMat = bst_history('add', DataMat, 'moveraw', ['Moved from ' oldPath ' to ' newPath]);
            % Save new file structure
            DataMat.F    = sFile;
            bst_save(DataFile, DataMat, 'v6');
        end
        
        OutputFiles{end+1} = sInputs(i).FileName;
    end
end

    
%% ===== COMPUTE =====
function [sFile, warnMessage] = Compute(sFile, oldPath, newPath)
    % ===== PARSE INPUTS =====
    if ~isstruct(sFile)
        DataFile = file_fullpath(sFile);
        DataMat = in_bst_data(DataFile);
        sFile = DataMat.F;
    end
    warnMessage = [];
    
    % File name
    [fPath, fBase, fExt] = bst_fileparts(sFile.filename);
    BaseName = [fBase, fExt];
    
    % Replace for new path if applicable
    if strfind(fPath, oldPath)
        fNewPath = strrep(fPath, oldPath, newPath);
        newFilename = bst_fullfile(fNewPath, BaseName);
        
        if exist(newFilename, 'file') == 2
            sFile.filename = newFilename;
            
            % Update the list of MEG4 files, in case of CTF recordings            
            if strcmpi(sFile.format, 'CTF') || strcmpi(sFile.format, 'CTF-CONTINUOUS')
                for i = 1:length(sFile.header.meg4_files)
                    [tmp, fBase, fExt] = bst_fileparts(sFile.header.meg4_files{i});
                    sFile.header.meg4_files{i} = bst_fullfile(fNewPath, [fBase, fExt]);
                end
            end
        else
            warnMessage = ['Invalid new path "' fNewPath '" for file "' sFile.filename '".'];
        end
    else
        warnMessage = ['Pattern not found in file "' sFile.filename '".'];
    end    
end


