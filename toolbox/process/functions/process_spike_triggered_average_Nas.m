function varargout = process_spike_triggered_average( varargin )
% PROCESS_SPIKE_TRIGGERED_AVERAGE: Computes the spike triggered average.
% 

% There are two different TimeWindow Notations here:
% 1. Timewindow around the spike (This is the one that is asked as input when the function is called).
% 2. Timewindow of the trials imported to the function.

% The function selects a TimeWindow around the Spike of a specific neuron.
% Then averages the LFPs oe each electrode.
% If this Spike TimeWindow is outside the TimeWindow of the Trial, the
% spike is ignored for computation.



% USAGE:    sProcess = process_spike_field_coherence('GetDescription')
%        OutputFiles = process_spike_field_coherence('Run', sProcess, sInput)

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
% Authors: Konstantinos Nasiotis, 2018

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Spike Triggered Average';
    sProcess.FileTag     = 'STA';
    sProcess.Category    = 'custom';
    sProcess.SubGroup    = 'e-Phys Functions';
    sProcess.Index       = 2507;
    sProcess.Description = 'http://www.in.gr';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Options: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'EEG';   % MAYBE ADD RAW SUPPORT HERE AS WELL
    sProcess.options.sensortypes.InputTypes = {'data'};
    sProcess.options.sensortypes.Group   = 'input';
    % Options: Parallel Processing
    sProcess.options.paral.Comment = 'Parallel processing';
    sProcess.options.paral.Type    = 'checkbox';
    sProcess.options.paral.Value   = 1;
    % Options: Segment around spike
    sProcess.options.timewindow.Comment  = 'Spike Time window: ';
    sProcess.options.timewindow.Type     = 'range';
    sProcess.options.timewindow.Value    = {[-0.150, 0.150],'ms',[]};
   
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Initialize returned values
    OutputFiles = {};
    % Extract method name from the process name
    strProcess = strrep(strrep(func2str(sProcess.Function), 'process_', ''), 'timefreq', 'morlet');
    
    % Add other options
    tfOPTIONS.Method = strProcess;
    if isfield(sProcess.options, 'sensortypes')
        tfOPTIONS.SensorTypes = sProcess.options.sensortypes.Value;
    else
        tfOPTIONS.SensorTypes = [];
    end    
    
    % If a time window was specified
    if isfield(sProcess.options, 'timewindow') && ~isempty(sProcess.options.timewindow) && ~isempty(sProcess.options.timewindow.Value) && iscell(sProcess.options.timewindow.Value)
        tfOPTIONS.TimeWindow = sProcess.options.timewindow.Value{1};
    elseif ~isfield(tfOPTIONS, 'TimeWindow')
        tfOPTIONS.TimeWindow = [];
    end
    
    
    %%%%%%%%%%%%%%%%%% MARTIN - ARE THESE OUTPUTS CORRECT? %%%%%%%%%%%%%%%%
    % Output
    if isfield(sProcess.options, 'avgoutput') && ~isempty(sProcess.options.avgoutput) && ~isempty(sProcess.options.avgoutput.Value)
        if sProcess.options.avgoutput.Value
            tfOPTIONS.Output = 'average';
        else
            tfOPTIONS.Output = 'all';
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tfOPTIONS.TimeVector = in_bst(sInputs(1).FileName, 'Time');

    if sProcess.options.timewindow.Value{1}(1)>=0 || sProcess.options.timewindow.Value{1}(2)<=0
        bst_report('Error', sProcess, sInputs, 'The time-selection must be around the spikes.');
    elseif sProcess.options.timewindow.Value{1}(1)==tfOPTIONS.TimeVector(1) && sProcess.options.timewindow.Value{1}(2)==tfOPTIONS.TimeVector(end)
        bst_report('Error', sProcess, sInputs, 'The spike window has to be smaller than the trial window');
    end
 
    
    % === OUTPUT STUDY ===
    % Get output study
    [~, iStudy, ~] = bst_process('GetOutputStudy', sProcess, sInputs);
    tfOPTIONS.iTargetStudy = iStudy;
    
    % Get channel file
    sChannel = bst_get('ChannelForStudy', iStudy);
    % Load channel file
    ChannelMat = in_bst_channel(sChannel.FileName);
    
    
    % === START COMPUTATION ===
    sampling_rate = round(abs(1. / (tfOPTIONS.TimeVector(2) - tfOPTIONS.TimeVector(1))));
    
    [temp, ~] = in_bst(sInputs(1).FileName);
    
    nElectrodes = 0;
    for iChannel = 1:length(ChannelMat.Channel)
       if ChannelMat.Channel(iChannel).Type == 'EEG' % Maybe we can add this option to be available on the raw file as well???
          nElectrodes = nElectrodes + 1;               
       end
    end

    nTrials = length(sInputs);
    time_segmentAroundSpikes = linspace(sProcess.options.timewindow.Value{1}(1), sProcess.options.timewindow.Value{1}(2), abs(sProcess.options.timewindow.Value{1}(2))* sampling_rate + abs(sProcess.options.timewindow.Value{1}(1))* sampling_rate + 1);    

    
    % Prepare parallel pool, if requested
    if sProcess.options.paral.Value
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            parpool;
        end
    else
        poolobj = [];
    end
    
    
    
    
    %% Collect all the average LFPs for each trial for all Neurons.
    everything = struct(); % This is a struct 1xnTrials
    
    % I get the files outside of the parfor so it won't fail.
    % This loads the information from ALL TRIALS on ALL_TRIALS_files
    % (Shouldn't create a memory problem).
    ALL_TRIALS_files = struct();
    for iFile = 1:nTrials
        ALL_TRIALS_files(iFile).a = in_bst(sInputs(iFile).FileName);
    end
    
    
    % Optimize this
    if ~isempty(poolobj) 
        parfor iFile = 1:nTrials
            [LFPs_single_trial] = get_LFPs(ALL_TRIALS_files(iFile).a, nElectrodes, sProcess, time_segmentAroundSpikes, sampling_rate);
            everything(iFile).LFPs_single_trial = LFPs_single_trial;
        end 
    else
        for iFile = 1:nTrials
            [trial, ~] = in_bst(sInputs(iFile).FileName);
            [LFPs_single_trial] = get_LFPs(trial, nElectrodes, sProcess, time_segmentAroundSpikes, sampling_rate);
            everything(iFile).LFPs_single_trial = LFPs_single_trial;
        end 
    end
        
        
    
    %% Calculate the SFC
    % The Spike Field Coherence should be a 3d matrix
    % Number of neurons x Frequencies x Electrodes
    % Ultimately the user will select the NEURON that wants to be displayed,
    % and a 2D image with the other two dimensions will appear, showing the
    % coherence of the spikes of that neuron with the LFPs on every
    % electrode on all frequencies.
    

    % Create a cell that holds all of the labels and one for the unique labels
    % This will be used to take the averages using the appropriate indices
    all_labels = struct;
    labelsForDropDownMenu = {}; % Unique neuron labels (each trial might have different number of neurons). We need everything that appears.
    for iFile = 1:nTrials
        for iNeuron = 1:length(everything(iFile).LFPs_single_trial)
            all_labels.labels{iNeuron,iFile} = everything(iFile).LFPs_single_trial(iNeuron).label;
            labelsForDropDownMenu{end+1} = everything(iFile).LFPs_single_trial(iNeuron).label;
        end
    end
    all_labels = all_labels.labels;
    labelsForDropDownMenu = unique(labelsForDropDownMenu,'stable');
    
    
    
    
    
    
    
    STA = zeros(length(labelsForDropDownMenu), length(time_segmentAroundSpikes), nElectrodes); % Number of neurons x Frequencies x Electrodes : 161x301x192
    tempAverageLFP = zeros(nElectrodes, length(time_segmentAroundSpikes));                  % 192 x 301
    
    
    for iNeuron = 1:length(labelsForDropDownMenu)
        %% For each TRIAL, get the index of the label that corresponds to the appropriate neuron.
        logicalEvents = ismember(all_labels, labelsForDropDownMenu{iNeuron}); % Find the index of the spike-events that correspond to that electrode (Exact string match). This linearizes the cell. I need to dilenearize it.
        iEvents = zeros(size(all_labels,2),1);
        for iTrial = 1:size(all_labels,2)
            [temp, ~] = find(logicalEvents(:,iTrial));
            if ~isempty(temp)
                iEvents(iTrial) = temp;
            else
                iEvents(iTrial) = 0; % This shows that that neuron didn't fire any spikes on that trial
            end
        end
        
        
        %% Take the Averagesof the appropriate indices
        for iTrial = 1:size(all_labels,2)
        
            divideBy = 0;
            if iEvents(iTrial)~=0
                tempAverageLFP = tempAverageLFP + everything(iTrial).LFPs_single_trial(iEvents(iTrial)).nSpikes * everything(iTrial).LFPs_single_trial(iEvents(iTrial)).avgLFP; % The avgLFP are sum actually. 
                divideBy = divideBy + everything(iTrial).LFPs_single_trial(iEvents(iTrial)).nSpikes;
            end 
        end
        
        STA(iNeuron,:,:) = (tempAverageLFP./divideBy)'; % 
           
    end
    
%     
%     
%     %% Plot an example for proof of concept
%     iNeuron = 100;
%     figure(1);
%     imagesc(squeeze(SFC(iNeuron,:,:))')        % SFC: Number of neurons x Frequencies x Electrodes : 161x151x192
%     ylabel 'iElectrode'
%     xlabel 'Frequency (Hz)'
%     title ({'Spike Field Coherence';['Neuron ' num2str(iNeuron)]})
%     
%     
    


   
    %% Rename the labels to something more meaningfull
    labelsForDropDownMenu_temp = cell(1,length(labelsForDropDownMenu));
    for iNeuron = 1:length(labelsForDropDownMenu)
        multiple_neurons_indicator = strfind(labelsForDropDownMenu{iNeuron},'|');
        if isempty(multiple_neurons_indicator) % Only a single neuron on the electrode
            wordsInLabel = strsplit(labelsForDropDownMenu{iNeuron}, ' ');
            labelsForDropDownMenu_temp{iNeuron} = ['Electrode ' num2str(wordsInLabel{3}) ': Neuron 1']; % Rename 'Spikes Electrode 1' to Electrode 1 : Neuron 1
        else
            wordsInLabel = strsplit(labelsForDropDownMenu{iNeuron}, ' ');
            labelsForDropDownMenu_temp{iNeuron} = ['Electrode ' num2str(wordsInLabel{3}) ': Neuron ' num2str(labelsForDropDownMenu{iNeuron}(multiple_neurons_indicator(1)+1:multiple_neurons_indicator(2)-1))];
            % Rename 'Spikes Electrode 1 |2|' to Electrode 1 : Neuron 2
        end
    end
    labelsForDropDownMenu = labelsForDropDownMenu_temp;



    
    %%
    tfOPTIONS.ParentFiles = {sInputs.FileName};
    
    % Prepare output file structure
    FileMat.TF = STA;
    FileMat.Time = time_segmentAroundSpikes; % These values are in order to trick Brainstorm with the correct values (This needs to be improved. Talk to Martin)
    FileMat.TFmask = [];
    FileMat.Freqs = 1:nElectrodes;      % These values are in order to trick Brainstorm with the correct values (This needs to be improved. Talk to Martin)
    FileMat.Std = [];
    FileMat.Comment = 'Spike Triggered Average';
    FileMat.DataType = 'data';
    FileMat.TimeBands = [];
    FileMat.RefRowNames = [];
    FileMat.RowNames = labelsForDropDownMenu;
    FileMat.Measure = 'power';
    FileMat.Method = 'morlet';
    FileMat.DataFile = []; % Leave blank because multiple parents
    FileMat.SurfaceFile = [];
    FileMat.GridLoc = [];
    FileMat.GridAtlas = [];
    FileMat.Atlas = [];
    FileMat.HeadModelFile = [];
    FileMat.HeadModelType = [];
    FileMat.nAvg = [];
    FileMat.ColormapType = [];
    FileMat.DisplayUnits = [];
    FileMat.Options = tfOPTIONS;
    FileMat.History = [];
   
    
    % Get output study
    sTargetStudy = bst_get('Study', iStudy);
    % Output filename
    FileName = bst_process('GetNewFilename', bst_fileparts(sTargetStudy.FileName), 'timefreq_spike_field_coherence');
    OutputFiles = {FileName};
    % Save output file and add to database
    bst_save(FileName, FileMat, 'v6');
    db_add_data(tfOPTIONS.iTargetStudy, FileName, FileMat);
    % Display report to user
    bst_report('Info', sProcess, sInputs, 'Success');
    disp('BST> process_spike_field_coherence: Success');
    
    
    % Close parallel pool
    if sProcess.options.paral.Value
        if ~isempty(poolobj)
            delete(poolobj);
        end
    end
end






function all = get_LFPs(trial, nElectrodes, sProcess, time_segmentAroundSpikes, sampling_rate)


    allEventLabels = {trial.Events.label};
        
    %% Get the events that include neurons

    % Important Variable here!
    spikeEvents = []; % The spikeEvents variable holds the indices of the events that correspond to spikes.
    %%%

    for ielectrode = 1:nElectrodes

        % Check for single neuron on electrode
        iEvent = find(ismember(allEventLabels, ['Spikes Electrode ' num2str(ielectrode)])); % Find the index of the spike-events that correspond to that electrode (Exact string match)

        % Check for multiple neurons on the same electrode
        if isempty(iEvent)
            Index = find(contains(allEventLabels,['Spikes Electrode ' num2str(ielectrode)])); % Check that the spike string is contained (this will be true in electrodes with multiple neurons: Spikes Electrode 150 |1|, Spikes Electrode 150 |2|)
            for iNeuron = 1:length(Index)
                if strcmp(allEventLabels(Index(iNeuron)), ['Spikes Electrode ' num2str(ielectrode) ' |' num2str(iNeuron) '|'])
                    iEvent = Index(iNeuron);
                    spikeEvents(end+1) = iEvent;
                end 
            end
        else
            spikeEvents(end+1) = iEvent;
        end

    end
    clear iEvent ielectrode allEventLabels


    all = struct();
    %% Get segments around each spike, FOR EACH NEURON
    for iNeuron = 1:length(spikeEvents) % iNeuron is the iEvent

        % Check that the entire segment around the spikes [-150,150]ms
        % is inside the trial segment and keep only those events
        events_within_segment = trial.Events(spikeEvents(iNeuron)).samples(trial.Events(spikeEvents(iNeuron)).times > trial.Time(1)   + abs(sProcess.options.timewindow.Value{1}(1)) & ...
                                                                           trial.Events(spikeEvents(iNeuron)).times < trial.Time(end) - abs(sProcess.options.timewindow.Value{1}(2)));

        %% Create a matrix that holds all the segments around the spike
        % of that neuron, for all electrodes.
        allSpikeSegments_singleNeuron_singleTrial = zeros(length(events_within_segment),size(trial.F,1),abs(sProcess.options.timewindow.Value{1}(2))* sampling_rate + abs(sProcess.options.timewindow.Value{1}(1))* sampling_rate + 1);

        for ispike = 1:length(events_within_segment)
            allSpikeSegments_singleNeuron_singleTrial(ispike,:,:) = trial.F(:, round(length(trial.Time) / 2) + events_within_segment(ispike) - abs(sProcess.options.timewindow.Value{1}(1)) * sampling_rate: ...
                                                                               round(length(trial.Time) / 2) + events_within_segment(ispike) + abs(sProcess.options.timewindow.Value{1}(2)) * sampling_rate  ...
                                                                           );
        end

        all(iNeuron).label   = trial.Events(spikeEvents(iNeuron)).label;
        all(iNeuron).nSpikes = length(events_within_segment);
        all(iNeuron).avgLFP  = squeeze(sum(allSpikeSegments_singleNeuron_singleTrial,1));
        all(iNeuron).Used    = 0; % This indicates if this entry has already been used for computing the SFC (some spikes might not appear on every trial imported, so a new Neuron should be identified on a later trial).


    end
        
end