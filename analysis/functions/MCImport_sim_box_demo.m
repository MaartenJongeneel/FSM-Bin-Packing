function mc_rtc_data = MCImport_sim_box_demo(MCdatalog,meta_sens,timestamp);
% This function creates the mc_rtc_data folder for the SENSOR_MEASUREMENT
% folder. It contains the raw data logged by mc_rtc and the
% POSTPROCESSSING folder containing the postprocessed data.
%
% INPUTS:    MCdatalog   : .csv file of a recording (mc_rtc datalog)
%            SETUPINFO   : path to folder containing all information of
%                          the objects, environment, and robots. In the
%                          SETUPINFO folder, three folders should be
%                          present: OBJECT, ROBOT and ENVIRONMENT. 
%            meta        : struct containing metadata of the recording
%
% OUTPUTS:   mc_rtc_data : struct containing the raw data and
%                          POSTPROCESSING data of the mc_rtc data
%% Load mc_rtc data and write it to a struct
%%
%Get the creation timestamp from the file name
% timestamp = extractAfter(extractBefore(MCdatalog,'.csv'),'mc_log_'); %Get the time of creation

if nargin == 1
    meta_sens = [];
    timestamp = [];
end
%% Setup the main data-import
%Set the Text Import Options
opts = detectImportOptions(MCdatalog);
opts.SelectedVariableNames = {'t'};
idx = 2;
for ii = 1:size(opts.VariableNames,2)
    %select the position and orientation of all the boxes
    if contains(opts.VariableNames{ii},'box',IgnoreCase=true) && (contains(opts.VariableNames{ii},'FloatingBase_position',IgnoreCase=true) || contains(opts.VariableNames{ii},'FloatingBase_orientation',IgnoreCase=true))
        opts.SelectedVariableNames{idx} = opts.VariableNames{ii};
        idx = idx + 1;
%     else
%         if ~contains(opts.VariableNames{ii},'box',IgnoreCase=true) && ~contains(opts.VariableNames{ii},'UR10FTSensor',IgnoreCase=true) && ~contains(opts.VariableNames{ii},'tote1',IgnoreCase=true) && ~contains(opts.VariableNames{ii},'tote2',IgnoreCase=true)
%             fprintf(opts.VariableNames{ii})
%             fprintf('\n')
%         end 
    end
    %select the wrench measuremets from the Force/Torque sensor
    if contains(opts.VariableNames{ii},'UR10FTSensor',IgnoreCase=true) 
        opts.SelectedVariableNames{idx} = opts.VariableNames{ii};
        idx = idx + 1;
    end
    
    if contains(opts.VariableNames{ii},'tote1',IgnoreCase=true) && (contains(opts.VariableNames{ii},'position',IgnoreCase=true) || contains(opts.VariableNames{ii},'orientation',IgnoreCase=true))
        opts.SelectedVariableNames{idx} = opts.VariableNames{ii};
        idx = idx + 1;
    end
end

%Import the data
ds = readtable(MCdatalog,opts);

%Change time indication Header to "Time(s)"
ds.Properties.VariableNames{string(ds.Properties.VariableNames)' == 't'} = 'Time(s)';

%% POSTPROCESSING
ii = 1;
while(1)
    % get the time
    if strcmp(opts.SelectedVariableNames{ii},'t') 
        mc_rtc_data.time = table2array(ds(:,1));
        ii = ii + 1;
    end
    % postprocess FT sensor data
    if contains(opts.SelectedVariableNames{ii},'UR10FTSensor_c',IgnoreCase=true) 
        F = table2array(ds(:,ii+3:ii+5));
        tau = table2array(ds(:,ii:ii+2));
        mc_rtc_data.FT_sensor = [F';tau'];
        ii = ii + 6;
    end
    % postprocess items pose data
    if contains(opts.SelectedVariableNames{ii},'box',IgnoreCase=true) || contains(opts.SelectedVariableNames{ii},'tote1',IgnoreCase=true)
        itemName = extractBefore(opts.SelectedVariableNames{ii},'_FloatingBase');
        if isempty(itemName)
            itemName = 'tote1'; % if is not a box item, is the tote.. extend with other items if you want to include them here
        end
        q = table2array(ds(:,ii:ii+3));
        p = table2array(ds(:,ii+4:ii+6));
        H = [];
        for jj = 1:size(q,1)
            H(:,:,jj) = [quat2rotm(q(jj,:)) p(jj,:)'; 0 0 0 1];
        end
        mc_rtc_data.(itemName).transform = H;
        ii = ii + 7;
    end
    
    if ii >= size(opts.SelectedVariableNames,2)
        break
    end

end


%% POSTPROCESSING
%Get the transformation matrices describing the pose of each link w.r.t.
%the base. This procedured is based on the the URDF of the robot.
% geompath = fullfile('writeH5/setupinfo/ROBOT/Panda002/Description/');
% data.URDF.ds = fileread(strcat(geompath,'Panda002.urdf'));

% geompath = fullfile('./writeH5/setupinfo/ROBOT/Manipulator/UR10/Description/');
% data.URDF.ds = fileread(strcat(geompath,'robot.txt'));

%Get name per link from urdf file
% linknames = string(regexp(data.URDF.ds,'<link name="(.*)"','tokens','dotexceptnewline'));

%Load the urdf file to create a rigidbodytree
% olddir = cd(geompath);
% bot = importrobot(data.URDF.ds);
% cd(olddir)

%Set default base position to origin without rotations (otherwise base is 1m high)
% setFixedTransform(bot.Bodies{1}.Joint,eye(4));

%Change dataformat to use rows
% bot.DataFormat = 'row';

%Obtain the actual joint position data from mc_rtc_data
% panda1_jointpos = startsWith(string(mc_rtc_data.datalog.ds.Properties.VariableNames)','panda_1_q_');
% panda2_jointpos = startsWith(string(mc_rtc_data.datalog.ds.Properties.VariableNames)','panda_2_q_');
% q_panda1 = table2array(mc_rtc_data.datalog.ds(:,panda1_jointpos));
% q_panda2 = table2array(mc_rtc_data.datalog.ds(:,panda2_jointpos));

%Calculation of transform per body per time step
% for ii = 1:ceil(height(mc_rtc_data.datalog.ds))
%     for il = 1:length(linknames(2:end)) %for each link
%         if contains(lower(linknames(il)),'panda_link0')
%             H1 = eye(4); 
% %             H2 = eye(4);
%         else
%             try
%                 H1 = getTransform(bot,q_panda1(ii,:),bot.BodyNames{il},'panda_link0');
%             catch
%                 continue;
%             end
% %             H2 = getTransform(bot,q_panda2(ii,:),bot.BodyNames{il},'panda_link0');
%         end
%         mc_rtc_data.POSTPROCESSING.FrankaEmikaRobot.(linknames(il)).transforms.ds{ii} = H1;
% %         mc_rtc_data.POSTPROCESSING.Panda003.(linknames(il)).transforms.ds{ii} = H2;
%     end    
%     %textwaitbar(ii,ceil(height(mc_rtc_data.datalog.ds)),'Importing data mc_rtc      ');
% end

%Add metadata to the transformation matrices per link
% for il = 1:length(linknames)
%     mc_rtc_data.POSTPROCESSING.FrankaEmikaRobot.(linknames(il)).transforms.attr.note = strcat("4x4 transformation matrices over time of ", linknames(il), " expressing its pose in terms of the base of the robot.");
% %     mc_rtc_data.POSTPROCESSING.Panda003.(linknames(il)).transforms.attr.note = strcat("4x4 transformation matrices over time of ", linknames(il), " expressing its pose in terms of the base of the robot.");
% end

%% Add metadata
% mc_rtc_data.datalog.attr = meta_sens.datalog.attr; 
% mc_rtc_data.datalog.attr.duration = {num2str(table2array(mc_rtc_data.datalog.ds(end,1)))};
% mc_rtc_data.datalog.attr.duration_unit = {'s'};
% mc_rtc_data.datalog.attr.creation_timestamp = {char(timestamp)}; 
% mc_rtc_data.datalog.attr.sample_frequency = {num2str((height(mc_rtc_data.datalog.ds)-1)/table2array(mc_rtc_data.datalog.ds(end,1)))};
% mc_rtc_data.datalog.attr.sample_frequency_unit = {'Hz'};
% mc_rtc_data.attr.model={'Franka Emika Panda embedded sensors'};
% mc_rtc_data.attr = meta_sens.attr;
% mc_rtc_data.POSTPROCESSING.attr.note = {'Folder containing transformation matrices of the different parts of the robot, computed from datalog'};

%% Final checks
% nanvals = any(isnan(table2array(mc_rtc_data.datalog.ds)),2);
% if any(nanvals)
%     warning(["MCImport: "+ sum(nanvals)+" bad measurements in mc_rtc_data.datalog.ds! ("+100*sum(nanvals)/mc_rtc_data.samples+"%)"]);
% end

end
