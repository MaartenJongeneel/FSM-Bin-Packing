%% Genaral Configuration
set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
close all; clearvars; clc;
addpath( genpath( pwd ) );

%% Settings
seq = 3;     %Select the sequence you want to analyze 
verbose = true;   %If you want to print figures and text 

% fig_xsize = 600;
% doSave = 0;

%% 1)Process collected raw experimental data
ArchiveFileName = convertStringsToChars(append('data/Box-demonstratorSeq',string(seq),'.h5'));         %Select the experimental data
fn_sim = append('data/Seq',string(seq),'/sim_data.csv'); %Select the simulation data

% Load the H5 Archive
fprintf( 'Loading raw data from %s \nThis may take some time... \n', ArchiveFileName);
Data = readH5(ArchiveFileName);

% Retrive and plot some information about the data in the Archive
fields = fieldnames(Data);
data_length = length(fields);
if verbose
    fprintf( [' The Archive: \"' ArchiveFileName '\" contains ' num2str(data_length-2) ' experiments. \n\n'] );
    for idx=3:data_length
        if isstruct(Data.(fields{idx}))
            
            CStr =  cellstr(['Exp ', num2str(idx-2,'%03d'), ': ' ,fields{idx}, ' -> ', Data.(fields{idx}).attr.note] );
            fprintf( ' %s', CStr{:});
            fprintf('\n')
        end
    end
end 

%% 2) Process FT calibration dataset

ft_wrench = Data.(fields{3}).SENSOR_MEASUREMENT.ft_data.datalog.wrench.ds;
ft_imu = Data.(fields{3}).SENSOR_MEASUREMENT.ft_data.datalog.imu.ds;
ft_time = Data.(fields{3}).SENSOR_MEASUREMENT.ft_data.datalog.time.ds;

c_idx = find(abs(ft_imu(1,:)) > 0.15 & abs(ft_imu(2,:)) > 0.15 & abs(ft_imu(3,:)) < 9.7 , 1);

ft_time = ft_time - ft_time(c_idx);
ft_time(1:c_idx-1) = [];
ft_imu(:,1:c_idx-1) = [];
ft_wrench(:,1:c_idx-1) = [];

% given the calibration sequence recorded in the UR10's teach pendant, the
% time intervals in which data must considered are the following
if seq ==1
    idxs = [2.5 4; 6.5 8; 10.5 12; 15 16.5; 25.8 27.3; 32 33.5 ];
elseif seq==2
    idxs = [2.5 4; 6.5 8; 10.5 12; 15 16.5; 25.8 27.3; 32 33.1 ];
elseif seq==3
    idxs = [2.5 4; 6.5 8; 10.5 12; 15 16.5; 25.8 27.3; 32 32.9 ];
elseif seq==4
    idxs = [2.5 4; 6.5 8; 10.5 12; 15 16.5; 25.8 27.3; 32 32.9 ];
elseif seq==5
    idxs = [2.5 4; 6.5 8; 10.5 12; 15 16.5; 25.8 27.3; 32 32.9 ];
end
% the calibration sequence has the main axis pointing downwards in the following order
% | 1 | 2 | 3 | 4 | 5 | 6 |
% |-x |+y |+x |-y |-z |+z |
bars = [];
for ii=1:size(idxs,1)
    init_idx = find(ft_time >= idxs(ii,1),1 );
    fin_idx  = find(ft_time >= idxs(ii,2),1 );
    bars = [bars ;init_idx; fin_idx]; 
    switch ii
        case 1
            f_(1,1) = mean(ft_wrench(1,init_idx:fin_idx));
            t_(1,4) = mean(ft_wrench(5,init_idx:fin_idx)); % t_y = -m.g.c_z 
            t_(2,6) = mean(ft_wrench(6,init_idx:fin_idx)); % t_z =  m.g.c_y
        case 2
            f_(2,2) = mean(ft_wrench(2,init_idx:fin_idx));
            t_(1,2) = mean(ft_wrench(4,init_idx:fin_idx)); % t_x = -m.g.c_z 
            t_(2,5) = mean(ft_wrench(6,init_idx:fin_idx)); % t_z =  m.g.c_x
        case 3
            f_(2,1) = mean(ft_wrench(1,init_idx:fin_idx));
            t_(2,4) = mean(ft_wrench(5,init_idx:fin_idx)); % t_y =  m.g.c_z 
            t_(1,6) = mean(ft_wrench(6,init_idx:fin_idx)); % t_z = -m.g.c_y
        case 4
            f_(1,2) = mean(ft_wrench(2,init_idx:fin_idx));
            t_(2,2) = mean(ft_wrench(4,init_idx:fin_idx)); % t_x =  m.g.c_z
            t_(1,5) = mean(ft_wrench(6,init_idx:fin_idx)); % t_z = -m.g.c_x
        case 5
            f_(1,3) = mean(ft_wrench(3,init_idx:fin_idx));
            t_(1,1) = mean(ft_wrench(4,init_idx:fin_idx)); % t_x = -m.g.c_y
            t_(2,3) = mean(ft_wrench(5,init_idx:fin_idx)); % t_y =  m.g.c_x 
            t_(1,7) = mean(ft_wrench(6,init_idx:fin_idx)); % t_z =  0
        case 6
            f_(2,3) = mean(ft_wrench(3,init_idx:fin_idx));
            t_(2,1) = mean(ft_wrench(4,init_idx:fin_idx)); % t_x =  m.g.c_y
            t_(1,3) = mean(ft_wrench(5,init_idx:fin_idx)); % t_y = -m.g.c_x 
            t_(2,7) = mean(ft_wrench(6,init_idx:fin_idx)); % t_z =  0
    otherwise
        disp('other value')
    end

end

f_offset = mean(f_);
f_ = f_ - f_offset;
s = mean(f_(2,:))./f_(2,:);
m = mean(f_(2,:))/9.81;

t_offset = mean(t_);
t_offset = [t_offset(2) t_offset(4) t_offset(7)];
t_ = t_ - [t_offset(1)*ones(2,2) t_offset(2)*ones(2,2) t_offset(3)*ones(2,3)];

c = [(t_(1,3)-t_(2,3))/2 (t_(1,1)-t_(2,1))/2 mean(t_(2,[2 4]))] ./(m*9.81);

if verbose
    figure()
    subplot(2,1,1)
    plot(ft_time,s(1)*(ft_wrench(1,:)-f_offset(1)),'r.',ft_time, s(2)*(ft_wrench(2,:)-f_offset(2)),'g.',ft_time,s(3)*(ft_wrench(3,:)-f_offset(3)),'b.')
    hold on
    % plot(ft_time,ft_imu(1,:),'r--',ft_time,ft_imu(2,:),'g--',ft_time,ft_imu(3,:),'b--')
    for ii=1:size(bars,1)
        plot(ft_time(bars(ii))*[1 1], 200*[-1 1], 'k-' )
    end
    grid on
    ylim([-20 20])

    subplot(2,1,2)
    plot(ft_time,ft_wrench(4,:)-t_offset(1),'r.',ft_time, ft_wrench(5,:)-t_offset(2),'g.',ft_time,ft_wrench(6,:)-t_offset(3),'b.')
    hold on
    % plot(ft_time,ft_imu(1,:),'r--',ft_time,ft_imu(2,:),'g--',ft_time,ft_imu(3,:),'b--')
    for ii=1:size(bars,1)
        plot(ft_time(bars(ii))*[1 1], 200*[-1 1], 'k-' )
    end
    grid on
    ylim([-3 3])
    
    [x,y,z] = sphere;
    x = x*f_(2,1);
    y = y*f_(2,2);
    z = z*f_(2,3);
    figure
    surf(x+f_offset(1),y+f_offset(2),z+f_offset(3))
    alpha 0.5
    hold on
    [x,y,z] = sphere;
    x = x*f_(2,1)*s(1);
    y = y*f_(2,2)*s(2);
    z = z*f_(2,3)*s(3);
    surf(x,y,z)
    alpha 0.5
    plot3(ft_wrench(1,:), ft_wrench(2,:),ft_wrench(3,:),'r.')
    plot3(s(1)*(ft_wrench(1,:)-f_offset(1)), s(2)*(ft_wrench(2,:)-f_offset(2)),s(3)*(ft_wrench(3,:)-f_offset(3)),'b.')
    grid on
    axis equal
end

% the estimated payload parameters of the gripper are
fprintf( ['Estimated parameters: \n-Gripper:\n  +mass: ' num2str(m) ' [kg]\n  +CoM: [' num2str(c) '] [m]\n'] );
fprintf( ['-FT sensor:\n  +force offset: [' num2str(f_offset) '] [N]\n  +torque offset: [' num2str(t_offset) '] [Nm]\n'] );

clear ft_wrench ft_time ft_imu
%% retrieve the data from the sequence
fprintf( ['Retreiving Mocap and F/T sensor data ...\n'] );
for ii = 4:6
    % Mocap data
    MocapFields = fieldnames(Data.(fields{ii}).SENSOR_MEASUREMENT.Mocap.POSTPROCESSING);
    Time{ii-3} = Data.(fields{ii}).SENSOR_MEASUREMENT.Mocap.datalog.ds{:,2};
    RBT_base{ii-3} = averageSE3(cat(3, Data.(fields{ii}).SENSOR_MEASUREMENT.Mocap.POSTPROCESSING.Base005.transforms.ds{:}));
    Tote{ii-3} = averageSE3(cat(3, Data.(fields{ii}).SENSOR_MEASUREMENT.Mocap.POSTPROCESSING.Tote001.transforms.ds{:}));
    Zivid{ii-3} = averageSE3(cat(3, Data.(fields{ii}).SENSOR_MEASUREMENT.Mocap.POSTPROCESSING.Zivid001.transforms.ds{:}));
    Fitting{ii-3} = ReconstructAndFilter( cat(3, Data.(fields{ii}).SENSOR_MEASUREMENT.Mocap.POSTPROCESSING.FittingCollar001.transforms.ds{:}),Time{ii-3},9, 20, 20, 2 );
    SC{ii-3} = cat(3, Data.(fields{ii}).SENSOR_MEASUREMENT.Mocap.POSTPROCESSING.SuctionCup004.transforms.ds{:}) ;
%     SC_ = ReconstructAndFilter( cat(3, Data.(fields{ii}).SENSOR_MEASUREMENT.Mocap.POSTPROCESSING.SuctionCup004.transforms.ds{:}) ,Time,100, 20, 20, 2 );
    for jj = 1:size(MocapFields,1)
        if startsWith(MocapFields{jj},'Box')
            Boxes.(MocapFields{jj}){ii-3} = cat(3, Data.(fields{ii}).SENSOR_MEASUREMENT.Mocap.POSTPROCESSING.(MocapFields{jj}).transforms.ds{:} );
        end
    end
    % FT sensor data
    % correct for the measurement offsets (use this for the other recordings in this sequence)
    ft_wrench{ii-3} = [diag(s) zeros(3,3); zeros(3,3), eye(3)]*(Data.(fields{ii}).SENSOR_MEASUREMENT.ft_data.datalog.wrench.ds - [f_offset';t_offset']);
    ft_imu{ii-3} = Data.(fields{ii}).SENSOR_MEASUREMENT.ft_data.datalog.imu.ds;
    ft_time{ii-3} = Data.(fields{ii}).SENSOR_MEASUREMENT.ft_data.datalog.time.ds;
    
    % UR10 data
    UR10.tcp{ii-3} = cat(3, Data.(fields{ii}).SENSOR_MEASUREMENT.UR10_sensor.POSTPROCESSING.TCP.transforms.ds{:});
    UR10.time{ii-3}=  Data.(fields{ii}).SENSOR_MEASUREMENT.UR10_sensor.datalog.ds.("Time(s)");
    
end

% express all the quantities in robot's base frame as is done for the data
% logged from mc_rtc in the simulation
% or convert the other way around, so from simulation data expressend in
% the robot's base frame into the motive world frame.
for jj = 1:3
    Tote{jj} =  RBT_base{jj}\Tote{jj};
    for ii=1:size(Time{jj},1)        
        if seq==1       
            Data_exp.box013{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box013{jj}(:,:,ii);
            Data_exp.box014{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box014{jj}(:,:,ii); 
            Data_exp.box015{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box015{jj}(:,:,ii); 
            Data_exp.box016{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box016{jj}(:,:,ii); 
            Data_exp.box017{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box017{jj}(:,:,ii); 
            Data_exp.box018{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box018{jj}(:,:,ii);
            Data_exp.box023{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box023{jj}(:,:,ii);
            Data_exp.box024{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box024{jj}(:,:,ii);
            Data_exp.box027{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box027{jj}(:,:,ii);
            Data_exp.box028{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box028{jj}(:,:,ii);
            Data_exp.box029{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box029{jj}(:,:,ii);
            Data_exp.box030{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box030{jj}(:,:,ii);
        elseif seq==2
            Data_exp.box014{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box014{jj}(:,:,ii); 
            Data_exp.box015{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box015{jj}(:,:,ii); 
            Data_exp.box017{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box017{jj}(:,:,ii); 
            Data_exp.box018{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box018{jj}(:,:,ii);
            Data_exp.box019{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box019{jj}(:,:,ii);
            Data_exp.box023{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box023{jj}(:,:,ii);
            Data_exp.box024{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box024{jj}(:,:,ii);
            Data_exp.box025{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box025{jj}(:,:,ii);
            Data_exp.box026{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box026{jj}(:,:,ii);
            Data_exp.box027{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box027{jj}(:,:,ii);
            Data_exp.box028{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box028{jj}(:,:,ii);
            Data_exp.box029{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box029{jj}(:,:,ii);
        elseif seq==3
            Data_exp.box013{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box013{jj}(:,:,ii); 
            Data_exp.box016{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box016{jj}(:,:,ii); 
            Data_exp.box017{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box017{jj}(:,:,ii); 
            Data_exp.box018{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box018{jj}(:,:,ii);
            Data_exp.box019{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box019{jj}(:,:,ii);
            Data_exp.box020{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box020{jj}(:,:,ii);
            Data_exp.box023{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box023{jj}(:,:,ii);
            Data_exp.box024{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box024{jj}(:,:,ii);
            Data_exp.box025{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box025{jj}(:,:,ii);
            Data_exp.box026{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box026{jj}(:,:,ii);
            Data_exp.box027{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box027{jj}(:,:,ii);
            Data_exp.box028{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box028{jj}(:,:,ii);
            Data_exp.box030{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box030{jj}(:,:,ii);
        elseif seq==4
            Data_exp.box017{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box017{jj}(:,:,ii);
            Data_exp.box023{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box023{jj}(:,:,ii); 
            Data_exp.box018{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box018{jj}(:,:,ii); 
            Data_exp.box030{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box030{jj}(:,:,ii); 
            Data_exp.box024{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box024{jj}(:,:,ii); 
            Data_exp.box025{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box025{jj}(:,:,ii);
            Data_exp.box019{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box019{jj}(:,:,ii);
            Data_exp.box014{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box014{jj}(:,:,ii);
            Data_exp.box026{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box026{jj}(:,:,ii);
            Data_exp.box020{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box020{jj}(:,:,ii);
            Data_exp.box029{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box029{jj}(:,:,ii);
            Data_exp.box015{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box015{jj}(:,:,ii); 
            Data_exp.box021{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box021{jj}(:,:,ii); 
            Data_exp.box022{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box022{jj}(:,:,ii); 
        elseif seq==5
            Data_exp.box013{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box013{jj}(:,:,ii);
            Data_exp.box016{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box016{jj}(:,:,ii); 
            Data_exp.box026{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box026{jj}(:,:,ii); 
            Data_exp.box014{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box014{jj}(:,:,ii); 
            Data_exp.box023{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box023{jj}(:,:,ii); 
            Data_exp.box024{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box024{jj}(:,:,ii);
            Data_exp.box025{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box025{jj}(:,:,ii);
            Data_exp.box017{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box017{jj}(:,:,ii);
            Data_exp.box015{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box015{jj}(:,:,ii);
            Data_exp.box018{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box018{jj}(:,:,ii);
            Data_exp.box027{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box027{jj}(:,:,ii);
            Data_exp.box019{jj}(:,:,ii) = RBT_base{jj}\Boxes.Box019{jj}(:,:,ii);
        end
    end
end
fprintf( ['Retreiving Mocap and F/T sensor data DONE\n'] );
%% Get simulation data
Data_sim = MCImport_sim_box_demo(fn_sim);
FT_sim = -Data_sim.FT_sensor;

%Correct for offset in z-direction for placement of Mocap frame w.r.t. simulation frame
for ii=1:length(Data_sim.time)
    Data_sim.box013.transform(:,:,ii) = Data_sim.box013.transform(:,:,ii)*[eye(3,3), [0; 0; -0.046]; zeros(1,3), 1]; 
    Data_sim.box014.transform(:,:,ii) = Data_sim.box014.transform(:,:,ii)*[eye(3,3), [0; 0; -0.046]; zeros(1,3), 1];
    Data_sim.box015.transform(:,:,ii) = Data_sim.box015.transform(:,:,ii)*[eye(3,3), [0; 0;  0.000]; zeros(1,3), 1];
    Data_sim.box016.transform(:,:,ii) = Data_sim.box016.transform(:,:,ii)*[eye(3,3), [0; 0; -0.010]; zeros(1,3), 1];
    Data_sim.box023.transform(:,:,ii) = Data_sim.box023.transform(:,:,ii)*[eye(3,3), [0; 0; -0.012]; zeros(1,3), 1];
end
  
%% Do a time alignment of Mocap and Simulation Data
%Here we choose to do a time alignment on Box018, which is present in each
%sequence. Since each sequence is recorded 3 times, you can choose which
%experiment you want to use for comparison. In the end, we will do the
%time-alignment based on the z-position of this box over time.
if seq == 1
    exp = 3;        %Select the experiment you want to compare with (1,2, or 3)
    t_vec = -0.05:0.01:0.05; %Time vector used to shift the simulation data
    comp = 2;       %Component of the vector for comparison
elseif seq == 2
    exp = 1;        %Select the experiment you want to compare with (1,2, or 3)
    t_vec = -0.6:0.05:-0.4; %Time vector used to shift the simulation data
    comp = 3;       %Component of the vector for comparison
elseif seq == 3
    exp = 3;        %Select the experiment you want to compare with (1,2, or 3)
    t_vec = -0.14:0.01:-0.10; %Time vector used to shift the simulation data
    comp = 3;       %Component of the vector for comparison
elseif seq == 4
    exp = 1;        %Select the experiment you want to compare with (1,2, or 3)
    t_vec = -0.6:0.05:-0.4; %Time vector used to shift the simulation data
    comp = 2;       %Component of the vector for comparison
elseif seq == 5
    exp = 3;        %Select the experiment you want to compare with (1,2, or 3)
    t_vec = -0.05:0.005:0.05; %Time vector used to shift the simulation data
    comp = 1;       %Component of the vector for comparison
end

fprintf( ['Performing a time alignment Mocap and Sim for Sequence ' num2str(seq) ' ...\n'] );

%Obtain the positions used for comparison
pos_moc = squeeze(Data_exp.box018{exp}(comp,4,:)); %Positon from experiments
pos_sim = squeeze(Data_sim.box018.transform(comp,4,:));  %Position from simulation

%Obtain the time vectors used for comparison
t_moc = Time{exp};       %Time of motion capture
t_sim = Data_sim.time;   %Time of simulation

[indx_moc2sim,indx_ts_moc]=TimeAlign(t_vec,t_moc,t_sim,pos_moc,pos_sim,125);

%Plot the result
figure; plot(t_moc,pos_moc); hold on;plot(t_moc,pos_sim(indx_moc2sim(indx_ts_moc,:)))
clear pos_sim
fprintf( ['Performing a time alignment Mocap and Sim: found at t=' num2str(t_vec(indx_ts_moc)) ' [s].\n'] );
%% Do a time alignment of FT data and Simulation Data
%Here we choose to do a time alignment on the x-component of the force
if seq == 1
    t_vec = -0.45:0.01:-0.35; %Time vector used to shift the simulation data
    comp = 1;       %Component of the vector for comparison
elseif seq == 2
    t_vec = -0.1:0.05:0.1; %Time vector used to shift the simulation data
    comp = 1;       %Component of the vector for comparson
elseif seq == 3
    t_vec = -0.5:0.01:-0.3; %Time vector used to shift the simulation data
    comp = 1;       %Component of the vector for comparson
elseif seq == 4
    t_vec = 0.09:0.005:0.11; %Time vector used to shift the simulation data
    comp = 1;       %Component of the vector for comparson
elseif seq == 5
    t_vec = -0.5:0.01:-0.3; %Time vector used to shift the simulation data
    comp = 1;       %Component of the vector for comparson
end

fprintf( ['Performing a time alignment FT data and Simulation ...\n'] );

%Obtain the positions used for comparison
ft_exp = ft_wrench{exp}; %F/T data from experiments
ft_sim = FT_sim;         %F/T data from simulation

%Obtain the time vectors used for comparison
t_fts = ft_time{exp}';   %Time of force/torque sensor
t_sim = Data_sim.time;  %Time of simulation

[indx_sim2fts,indx_ts_fts]=TimeAlign(t_vec,t_sim,t_fts,ft_sim(comp,:)',ft_exp(comp,:)',500);

%Plot the result
figure; plot(t_sim,ft_sim(1,:),'r.-'); hold on;plot(t_sim,ft_exp(1,indx_sim2fts(indx_ts_fts,:)),'b.-')

fprintf( ['Performing a time alignment FT data and Simulation: found at t=' num2str(t_vec(indx_ts_fts)) ' [s].\n'] );
%% Analysis of sim-to-real
% Sequence 1
if seq == 1
    %  Box013 | Box014 | Box015 | Box016 | Box017 | Box018 | Box023 | Box024 | Box027 | Box028 | Box029 | Box030 |
    %For the boxes in this sequence, define the time windows of interest (these are the moments the box is moving)
    t_int = [45, 63, 112, 129,  93, 177, 219, 236,  2, 23, 159, 201;
             57, 88, 125, 158, 106, 198, 230, 250, 20, 41, 175, 218];   
end
% Sequence 2
if seq == 2
    %  Box014 | Box015 | Box017 | Box018 | Box019 | Box023 | Box024 | Box025 | Box026 | Box027 | Box028 | Box029 |
    %For the boxes in this sequence, define the time windows of interest (these are the moments the box is moving)
    t_int = [160, 175, 127, 144, 195, 210, 3, 28, 53, 78, 93, 111;
             172, 190, 137, 156, 206, 230, 25, 48, 74, 87, 105, 120]; 
end
% Sequence 3
if seq == 3
    %  Box013 | Box016 | Box017 | Box018 | Box019 | Box020 | Box023 | Box024 | Box025 | Box026 | Box027 | Box028 | Box030
    %For the boxes in this sequence, define the time windows of interest (these are the moments the box is moving)
    t_int = [ 94, 208, 32, 53, 131, 256, 227,  3, 150, 179, 114, 272, 75;
             109, 224, 49, 71, 145, 268, 250, 28, 173, 202, 127, 289, 89];
end
% Sequence 4
if seq == 4
    %  Box017 | Box023 | Box018 | Box030 | Box024 | Box025 | Box019 | Box014 | Box026 | Box020 | Box029 | Box015 | Box021 | Box022 |
    %For the boxes in this sequence, define the time windows of interest (these are the moments the box is moving)
    t_int = [2, 24, 42, 62, 83, 101, 118, 137, 153, 182, 200, 220, 240, 258;
             18, 35, 56, 76, 93, 112, 130, 147, 178, 194, 214, 234, 250, 272];
end
% Sequence 5
if seq == 5
    %  Box013 | Box016 | Box026 | Box014 | Box023 | Box024 | Box025 | Box017 | Box015 | Box018 | Box027 | Box019 |
    %For the boxes in this sequence, define the time windows of interest (these are the moments the box is moving)
    t_int = [ 4, 22, 39, 68,  87, 118, 148, 177, 199, 218, 238, 271;
             15, 33, 63, 81, 112, 141, 173, 192, 209, 232, 247, 283];
end

%For each of these boxes, we need to post-process some data (unwrapping angles etc) to be able to compute rms values etc. We store them in a cell array
fn = fieldnames(Data_exp);
for ii = 1:length(fn)
    [pos_exp{ii},rot_exp{ii},pos_sim{ii},rot_sim{ii},err_pos{ii},err_rot{ii}]=ComparePoses(Data_exp.(fn{ii}){exp},t_moc,Data_sim.(fn{ii}).transform(:,:,indx_moc2sim(indx_ts_moc,:)),t_moc,t_int(:,ii),false,fn{ii});
    [err_f{ii}, err_t{ii}] = CompareWrenches(ft_exp(:,indx_sim2fts(indx_ts_fts,:)),t_sim,ft_sim,t_sim,t_int(:,ii));
    rms_pos(ii) = rms(err_pos{ii}(~isnan(err_pos{ii})));
    rms_rot(ii) = rms(err_rot{ii}(~isnan(err_rot{ii})));
end

%Display RMS errors
T = array2table([rms_pos*1000; rad2deg(rms_rot)],'VariableNames',fn','RowName',{'RMS pos error [mm]','RMS rot error [deg]'});
disp(T)


%Display rest-pose errors
tab_err = [];
for ii = 1:length(fn)
    tab_err = [tab_err, [err_pos{ii}(end)*1000;rad2deg(err_rot{ii}(end))]];
end
T = array2table(tab_err,'VariableNames',fn','RowName',{'Rest-pose pos error [mm]','Rest-pose rot error [deg]'});
disp(T)

%% Functions
function [pos_exp,rot_exp,pos_sim,rot_sim,err_pos,err_rot]=ComparePoses(H1,t1,H2,t2,xlims,doSave,fn_box)
%Quick function for plotting the poses.
%H1,t1 : Data from Mocap
%H2,t2 : Data from simulation (rotation axis are flipped)
%xlims : x-axis limits for region of interest
figure('rend','painters','Position',[400,260,560,600]);
ha = tight_subplot(4,2,[.07 .08],[.06 .01],[0.08 0.01]);  %[gap_h gap_w] [lower upper] [left right]
ax = ['x' ; 'y'; 'z' ];
% update = true;
RotVec{1}=NaN(3,length(H1));
RotVec{2}=NaN(3,length(H2));

for ii = 1:length(H1)
    PosVec{1}(:,ii) = H1(1:3,4,ii);
    RotVec{1}(:,ii) = vee(logm(H1(1:3,1:3,ii)));
    PosVec{2}(:,ii) = H2(1:3,4,ii);
    RotVec{2}(:,ii) = -vee(logm(H2(1:3,1:3,ii))); %flipping data because mc_rtc is weird..
end

%Find jumping moments of pi
for jj = 1:2
    jump = [0,vecnorm((diff(RotVec{jj}'))')>(2*pi-0.05)];
    jumped = false;
    for ii = 1:length(jump)
        if ~jumped && jump(ii) == 1
            jumped = true;
        elseif jumped && jump(ii) == 1
            jumped = false;
        end
        idx(ii) =  jumped;
    end
    
    %The moments the vectors switched around, we switch them back
    RotVec{jj}(:,idx) = -RotVec{jj}(:,idx);
end

%Obtain the data within the time limits
idx = t1>xlims(1)&t1<xlims(2);
pos_exp = PosVec{1}(:,idx);
pos_sim = PosVec{2}(:,idx);
rot_exp = RotVec{1}(:,idx);
rot_sim = RotVec{2}(:,idx);
t_vec = t1(idx);

%Compute the errors and remove outliers
err_pos = vecnorm(pos_exp-pos_sim);
% [~,idx] = rmoutliers(err_pos,"percentiles",[0 97.5]);
idx = err_pos>0.20;
err_pos(idx) = NaN;
% pos_exp(:,idx) = NaN;
err_rot = vecnorm(rot_exp-rot_sim);
% [~,idx] = rmoutliers(err_rot,"percentiles",[0 97.5]);
idx = err_rot>deg2rad(60);
err_rot(idx) = NaN;
% rot_exp(:,idx) = NaN;

%Plotting the data
for ii=0:2
%     subplot(4,2,ii*2+1)
    axes(ha(ii*2+1));
    plot(t_vec,pos_exp(ii+1,:),'.-','color','b'); hold on;
    plot(t_vec,pos_sim(ii+1,:),'.-','color','r');
    grid on
    xlabel('time [s]')
    ylabel(['$p_' ax(ii+1) '$ [m]'])
    xlim(xlims)

%     subplot(4,2,(ii+1)*2)
    axes(ha((ii+1)*2));
    plot(t_vec,rad2deg(rot_exp(ii+1,:)),'.-','color','b'); hold on;
    plot(t_vec,rad2deg(rot_sim(ii+1,:)),'.-','color','r');
    grid on
    xlabel('time [s]')
    ylabel(['$\theta_' ax(ii+1) '$ [deg]'])
    xlim(xlims)
end

%     subplot(4,2,7)
    axes(ha(7));
    plot(t_vec,err_pos,'.-','Color','#77AC30');
    hold on
    grid on
    xlabel('time [s]')
    ylabel('pos\_error [m]')
    xlim(xlims)

%     subplot(4,2,8)
    axes(ha(8));
    plot(t_vec,rad2deg(err_rot),'.-','Color','#77AC30');
    hold on
    grid on
    xlabel('time [s]')
    ylabel('rot\_error [deg]')
    xlim(xlims)

    if doSave; fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,append('Seq1_',fn_box,'_poses.pdf'),'-dpdf','-vector')
    end
end

function [f_err, t_err] =CompareWrenches(ft_real,t_real,ft_sim,t_sim,xlims)
ft_err = ft_real - ft_sim;
idx = (t_sim>=xlims(1)&t_sim<=xlims(2));
t = t_sim(idx);
f_err = vecnorm(ft_err(1:3,idx));
t_err = vecnorm(ft_err(4:6,idx));

figure('rend','painters','Position',[965,260,560,600]);
ha = tight_subplot(4,2,[.07 .08],[.06 .01],[0.08 0.01]);  %[gap_h gap_w] [lower upper] [left right]
axes(ha(1));
plot(t_real,ft_real(1,:),'b.');hold on;
plot(t_sim,ft_sim(1,:),'r.-');grid on;
ylabel('$f_x [N]$')
xlabel('Time [s]')
xlim(xlims)

axes(ha(3));
plot(t_real, ft_real(2,:),'b.');hold on;
plot(t_sim, ft_sim(2,:),'r.-');grid on;
ylabel('$f_y [N]$')
xlabel('Time [s]')
xlim(xlims)

axes(ha(5));
plot(t_real,ft_real(3,:),'b.');hold on;
plot(t_sim,ft_sim(3,:),'r.-');grid on;
ylabel('$f_z [N]$')
xlabel('Time [s]')
xlim(xlims)

axes(ha(7));
plot(t,f_err,'.-','color','#77AC30');grid on;
ylabel('$e_f [N]$')
xlabel('Time [s]')
xlim(xlims)

axes(ha(2));
plot(t_real,ft_real(4,:),'b.');hold on;
plot(t_sim,ft_sim(4,:),'r.-');grid on;
ylabel('$\tau_x [Nm]$')
xlabel('Time [s]')
xlim(xlims)

axes(ha(4));
plot(t_real, ft_real(5,:),'b.');hold on;
plot(t_sim, ft_sim(5,:),'r.-');grid on;
ylabel('$\tau_y [Nm]$')
xlabel('Time [s]')
xlim(xlims)

axes(ha(6));
plot(t_real,ft_real(6,:),'b.');hold on;
plot(t_sim,ft_sim(6,:),'r.-');grid on;
ylabel('$\tau_z [Nm]$')
xlabel('Time [s]')
xlim(xlims)

axes(ha(8));
plot(t,t_err,'.-','color','#77AC30');grid on;
ylabel('$e_{\tau} [Nm]$')
xlabel('Time [s]')
xlim(xlims)
end

function [indx_hf2lf,indx_opt]=TimeAlign(t_vec,t_lf,t_hf,data_lf,data_hf,Fhf)
%t_vec : Time vector over which we vary
%t_lf  : Low-Frequency time signal (nx1)
%t_hf  : High-Frequency time signal (mx1)
%data_lf : Low Frequency data signal (nx1)
%data_hf : High Frequency data signal (mx1)
%Fhf     : Frequency of High Frequency signal (1x1)
%
%indx_hf2lf : index showing which hf indexes are used for each lf index (for each element of t_vec)
%indx_opt : index of t_vec that gives the optimum value
cnt=1;
for ii = 1:length(t_vec)
    tel=1;
    for i = 1:length(t_lf)
        abs_dif = abs(t_lf(i) - (t_hf + t_vec(ii))); %Compute the time difference
        [~,indx_hf2lf(ii,tel)] = min(abs_dif); %Find the indices for which the time difference is minimal
        tel=tel+1;
        cnt=cnt+1;
        textwaitbar(cnt, length(t_lf)*length(t_vec),'Computing time alignment')
    end
    indx_con(:,ii) = abs(t_lf-(t_hf(indx_hf2lf(ii,:))+t_vec(ii))) < (1/Fhf)/2; %Do not consider indices for which the time difference is more than half the sampling time
    Nerror(ii) = norm(data_lf(indx_con(:,ii)) - data_hf(indx_hf2lf(ii,indx_con(:,ii)))); %For the indices we do consider, compute the error
end
[~,indx_opt] = min(Nerror); %Find the time index for which the error is minimal
end

%% Old code
% figure;
% subplot(3,2,1)
% plot(ft_time{2},ft_wrench{2}(1,:),'b.')
% hold on
% plot(T_sim,FT_sim(1,:),'r.-')
% grid on
% ylabel('$f_x [N]$')
% xlim([0 250])
% 
% subplot(3,2,3)
% plot(ft_time{2}, ft_wrench{2}(2,:),'b.')
% hold on
% plot(T_sim, FT_sim(2,:),'r.-')
% grid on
% ylabel('$f_y [N]$')
% xlim([0 250])
% 
% subplot(3,2,5)
% plot(ft_time{2},ft_wrench{2}(3,:),'b.')
% hold on
% plot(T_sim,FT_sim(3,:),'r.-')
% grid on
% ylabel('$f_z [N]$')
% xlabel('Time [s]')
% xlim([0 250])
% 
% subplot(3,2,2)
% plot(ft_time{2},ft_wrench{2}(4,:),'b.')
% hold on
% plot(T_sim,FT_sim(4,:),'r.-')
% grid on
% ylabel('$\tau_x [Nm]$')
% xlim([0 250])
% 
% subplot(3,2,4)
% plot(ft_time{2}, ft_wrench{2}(5,:),'b.')
% hold on
% plot(T_sim, FT_sim(5,:),'r.-')
% grid on
% ylabel('$\tau_y [Nm]$')
% xlim([0 250])
% 
% subplot(3,2,6)
% plot(ft_time{2},ft_wrench{2}(6,:),'b.')
% hold on
% plot(T_sim,FT_sim(6,:),'r.-')
% grid on
% ylabel('$\tau_z [Nm]$')
% xlim([0 250])

% Plot6DspaceTrajectory(Fitting{2},0.05,500)
% Plot6DtimeTrajectory(UR10.tcp{2},UR10.time{2}, '','b.')

% figure()
% PlotFrame(eye(4),'world',0.05)
% hold on
% PlotFrame(RBT_base{2},'RBT_{base}', 0.05)
% PlotFrame(Tote{2},'bin', 0.05)
% axis equal
% 
% %--------------------------------------------------------------------------
% figure;
% subplot(2,1,1)
% plot(ft_time{2},ft_wrench{2}(1,:),'r.',ft_time{2}, ft_wrench{2}(2,:),'g.',ft_time{2},ft_wrench{2}(3,:),'b.')
% grid on
% 
% subplot(2,1,2)
% plot(ft_time{2},ft_wrench{2}(4,:),'r.',ft_time{2}, ft_wrench{2}(5,:),'g.',ft_time{2},ft_wrench{2}(6,:),'b.')
% grid on


