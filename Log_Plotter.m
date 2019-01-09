function Log_Plotter ()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
clc
clear
close all
if isappdata(0,'file')
    rmappdata(0,'file');
end

% posible plots to pick from
% button varable name, button display name, Plot title, Plot Y axis label, Plot varable, Plot colums in varable, varable 2, colums in varable 2 ect
plots{1} = {'copter_att','Copter Attitude','Roll Angle - Desired vs actual','Angle (deg)','ATT',[3,4]};
plots{2} = {'copter_att','Copter Attitude','Pitch Angle - Desired vs actual','Angle (deg)','ATT',[5,6]};
plots{3} = {'copter_att','Copter Attitude','Yaw Angle - Desired vs actual','Angle (deg)','ATT',[7,8]};
plots{4} = {'copter_PID','Copter PID','Pitch PID','','PIDP',[4,5,6]};
plots{5} = {'copter_PID','Copter PID','Roll PID','','PIDR',[4,5,6]};
plots{6} = {'copter_PID','Copter PID','Yaw PID','','PIDY',[4,5,6]};
plots{7} = {'Airspeed','Airspeed','Airspeed','m/s','CTUN',10};
plots{8} = {'Baro','Barometic Alt','Barometic Alt','Altitude (m)','BARO',3,'BARO2',3};
plots{9} = {'Quadplane_alr','Q Altitude','Q Altitude','QTUN',[5,6]};
plots{10}= {'copter_rate','Copter Rates','Roll Rate - Desired vs actual','Angle (deg/s)','RATE',[3,4,5]};
plots{11}= {'copter_rate','Copter Rates','Pitch Rate - Desired vs actual','Angle (deg/s)','RATE',[6,7,8]};
plots{12}= {'copter_rate','Copter Rates','Yaw Rate - Desired vs actual','Angle (deg/s)','RATE',[9,10,11]};
plots{13}= {'RC_contorl_in','Control in','Control inputst','PWM','AETR',[3,4,5,6,7]};
plots{14}= {'Speed','Speed','Speed (m/s)','GPS',11,'CTUN',10};

% Defualts
file_name = 'Input File Name';
print_messages_default = 1;
save_param_file_default = 1;
open_param_file_default = 0;


%{


%plot_log('PWM Out','Angle (deg)',{RCOU,RCOU_label},[3,4],{RCIN,RCIN_label},[3,8])

%plot_log('Big Throttle','PWM',{RCIN,RCIN_label},[3]) % 3 = throttle % 5 = pitch % roll = 4 % yaw = 6
%plot_log('Big Pitch','PWM',{RCIN,RCIN_label},[5]) % 3 = throttle % 5 = pitch % roll = 4 % yaw = 6
%plot_log('V Tail PWM Out','Angle (deg)',{RCOU,RCOU_label},[6,7])
%plot_log('Ailerons PWM Out','Angle (deg)',{RCOU,RCOU_label},[3,10])

%plot_log('Flaps in','PWM',{RCIN,RCIN_label},[7]) % 3 = throttle % 5 = pitch % roll = 4 % yaw = 6
%plot_log('Flaps PWM Out','Angle (deg)',{RCOU,RCOU_label},[12,13])

%plot_log('Yaw in','PWM',{RCIN,RCIN_label},[6]) % 3 = throttle % 5 = pitch % roll = 4 % yaw = 6
%plot_log('Throttles PWM Out','Angle (deg)',{RCOU,RCOU_label},[8,9])
%plot_log('PWM Out','PWM',{RCOU,RCOU_label},[6,7,8,9,10])%,4,5,6,7])
%plot_log('PWM Out','PWM',{RCOU,RCOU_label},[6,9])%,4,5,6,7])
%plot_log('PWM Out','PWM',{RCOU,RCOU_label},[7,8])%,4,5,6,7])
%plot_log('PWM Out','PWM',{RCOU,RCOU_label},[10])%,4,5,6,7])

%plot_log('Control outputs','',{AETR,AETR_label},[3,4,5,6,7])
%plot_log('Aileron','',{AETR,AETR_label},[3])
%plot_log('Throttle','',{AETR,AETR_label},[5])
%plot_log('Pitch','',{AETR,AETR_label},[4])

%plot_log('GPS','Ground Speed (m/s)',{GPS,GPS_label},11,{CTUN,CTUN_label},10)
%plot_log('GPS','Ground Speed (m/s)',{GPS2,GPS2_label},11,{CTUN,CTUN_label},10)

%plot_log('Pitch PID','',{PIQP,PIQP_label},[4,5,6])
%plot(PIQP(:,2)/1e6,PIQP(:,4)+PIQP(:,5)+PIQP(:,6))

%plot_log('Roll PID','',{PIQR,PIQR_label},[4,5,6])
%plot(PIQR(:,2)/1e6,PIQR(:,4)+PIQR(:,5)+PIQR(:,6))

%plot_log('Yaw PID','',{PIQY,PIQY_label},[4,5,6])
%plot(PIQY(:,2)/1e6,PIQY(:,4)+PIQY(:,5)+PIQY(:,6))

%plot_log('Angle of Attack','AoA (deg)',{AOA,AOA_label},3)

%plot_log('Navigaiton',' ',{NTUN,NTUN_label},6)

%plot_log('Bat Current','Current (A)',{BAT,BAT_label},5)

%plot_log('Roll Rate','Rate',{RATE,RATE_label},[3,4])
%plot_log('Pitch Rate','Rate',{RATE,RATE_label},[6,7])

% HEX Motors
%%{
plot_log('HEX Motor PWM Out','PWM',{RCOU,RCOU_label},[3,4,5,6,7,8])%,4,5,6,7])

plot_log('HEX Motor PWM Out','PWM',{RCOU,RCOU_label},[11,12,13,14,15,16])%,4,5,6,7])
plot_log('Throttle in','PWM',{RCIN,RCIN_label},5)

%plot_log('ESC Log','RPM',{ESC1,ESC1_label},3,{ESC2,ESC2_label},3,{ESC3,ESC3_label},3,{ESC4,ESC4_label},3,{ESC5,ESC5_label},3,{ESC6,ESC6_label},3)
%plot_log('ESC Log','Voltage',{ESC1,ESC1_label},4,{ESC2,ESC2_label},4,{ESC3,ESC3_label},4,{ESC4,ESC4_label},4,{ESC5,ESC5_label},4,{ESC6,ESC6_label},4)
%plot_log('ESC Log','Current (A)',{ESC1,ESC1_label},5,{ESC2,ESC2_label},5,{ESC3,ESC3_label},5,{ESC4,ESC4_label},5,{ESC5,ESC5_label},5,{ESC6,ESC6_label},5)
%plot_log('ESC Log','TEMP (deg C)',{ESC1,ESC1_label},6,{ESC2,ESC2_label},6,{ESC3,ESC3_label},6,{ESC4,ESC4_label},6,{ESC5,ESC5_label},6,{ESC6,ESC6_label},6)

plot_log('Current Draw','Current (A)',{BAT,BAT_label},5)

plot_log('Roll Angle - Desired vs actual','Angle (deg)',{ATT,ATT_label},[3,4])
plot_log('Pitch Angle - Desired vs actual','Angle (deg)',{ATT,ATT_label},[5,6])
plot_log('Yaw Angle - Desired vs actual','Angle (deg)',{ATT,ATT_label},[7,8])

plot_log('Roll Rate - Desired vs actual','Angle (deg/s)',{RATE,RATE_label},[3,4])
plot_log('Pitch Rate - Desired vs actual','Angle (deg/s)',{RATE,RATE_label},[6,7])
plot_log('Yaw Rate - Desired vs actual','Angle (deg/s)',{RATE,RATE_label},[9,10])

%%}

% Bi Copter
%{
plot_log('Roll Angle - Desired vs actual','Angle (deg)',{ATT,ATT_label},[3,4])
plot_log('Pitch Angle - Desired vs actual','Angle (deg)',{ATT,ATT_label},[5,6])
plot_log('Yaw Angle - Desired vs actual','Angle (deg)',{ATT,ATT_label},[7,8])

plot_log('Roll Rate - Desired vs actual','Angle (deg/s)',{RATE,RATE_label},[3,4])
plot_log('Pitch Rate - Desired vs actual','Angle (deg/s)',{RATE,RATE_label},[6,7])
plot_log('Yaw Rate - Desired vs actual','Angle (deg/s)',{RATE,RATE_label},[9,10])

plot_log('ESC Log','RPM',{ESC1,ESC1_label},3,{ESC2,ESC2_label},3)
plot_log('ESC Log','Voltage',{ESC1,ESC1_label},4,{ESC2,ESC2_label},4)
plot_log('ESC Log','Current (A)',{ESC1,ESC1_label},5,{ESC2,ESC2_label},5)
plot_log('ESC Log','TEMP (deg C)',{ESC1,ESC1_label},6,{ESC2,ESC2_label},6)
%}
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs if you dont want to use GUI
% copy varable from defaults and over write
file_name = '10.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple(ish) GIU
load_in = @(src,event)Load_file(src,event,file_name);

handles.figure = figure;

handles.print_messages = uicontrol('Style',...
    'radiobutton',...
    'String','Print MSG''s',...
    'Position',[10 350 100 30],...
    'Value',print_messages_default,...
    'HandleVisibility','off');

handles.save_param_file = uicontrol('Style',...
    'radiobutton',...
    'String','Save Param File',...
    'Position',[100 350 100 30],...
    'Value',save_param_file_default,...
    'HandleVisibility','off');

handles.open_param_file = uicontrol('Style',...
    'radiobutton',...
    'String','Open Param File',...
    'Position',[200 350 100 30],...
    'Value',open_param_file_default,...
    'HandleVisibility','off');

handles.open_param_file_default = uicontrol('Style',...
    'pushbutton',...
    'String',file_name,...
    'Position',[10 300 200 30],...
    'HandleVisibility','off',...
    'Callback', @select_file);

handles.load_file = uicontrol('Style',...
    'pushbutton',...
    'String','Reload',...
    'Position',[210 300 100 30],...
    'HandleVisibility','off',...
    'Callback', load_in);

handles.mode_shadeing = uicontrol('Style',...
    'radiobutton',...
    'String','mode shading',...
    'Position',[10 250 100 30],...
    'Value',save_param_file_default,...
    'Enable','off',...
    'HandleVisibility','off');

handles.time_markers = uicontrol('Style',...
    'radiobutton',...
    'String','time markers',...
    'Position',[100 250 100 30],...
    'Value',save_param_file_default,...
    'Enable','off',...
    'HandleVisibility','off');

handles.time_crop = uicontrol('Style',...
    'radiobutton',...
    'String','time crop',...
    'Position',[200 250 100 30],...
    'Value',save_param_file_default,...
    'Enable','off',...
    'HandleVisibility','off');

plot_catagorys = cell(1,length(plots));
for n = 1:length(plots)
    plot_catagorys{n} = plots{n}{1,1};
end
[plot_catagorys, index] = unique(plot_catagorys,'stable');
Y_no = ceil((length(plot_catagorys) / 3)) - 1;
y_spacing = 25;
x_spacing = 100;
xloc = 0:x_spacing:x_spacing*2;
yloc = 0:-y_spacing:-y_spacing*Y_no;
[X,Y] = meshgrid(xloc,yloc);
X = X';
Y = Y';
X = X(:);
Y = Y(:);
for n = 1:length(index)
    handles.(plot_catagorys{n}) = uicontrol('Style',...
        'radiobutton',...
        'String',plots{index(n)}{2},...
        'Position',[X(n)+10 Y(n)+200 100 30],...
        'Enable','off',...
        'Value',1,...
        'HandleVisibility','off');
end

handles.plot = uicontrol('Style',...
    'pushbutton',...
    'String','Plot',...
    'Position',[200 10 100 30],...
    'Enable','off',...
    'HandleVisibility','off',...
    'Callback', @plot_figures);

handles.plot_close = uicontrol('Style',...
    'pushbutton',...
    'String','Close Plots',...
    'Position',[100 10 100 30],...
    'Enable','off',...
    'HandleVisibility','off',...
    'Callback', @plot_close);

setappdata(0,'handles', handles);
setappdata(0,'plots', plots);
setappdata(0,'plot_catagorys', plot_catagorys);



%{
return


global Mode_markes
global TimeMark
global colour
global MODELOOKUP
global Time_markers
global TimeCrop




TimeMark = [2213,2215];%(s) %Time Markers (vert lines)
%TimeCrop = [929,NaN]; % (s) Time to crop to
%TimeCrop = [2204,2217]; % (s) Time to crop to
%TimeCrop = [500,585];
TimeCrop = [140,165]; % (s) Time to crop to
TimeCrop = [0,nan]; % (s) Time to crop to

Time_markers = false;
Mode_markes = true;

% Param compare
if 0
    clearvars -except PARM
    PRAM2 = load('new','PARM');
    
    for n = 1:length(PARM)
        
        index = NaN;
        for m = 1:length(PRAM2.PARM)
            if  strcmp(PARM{n,1},PRAM2.PARM{m,1})
                index = m;
                break;
            end
        end
            
        % Ignore Parameters
        if strfind(PARM{n,1},'INS_')
        elseif strfind(PARM{n,1},'STAT_')
        elseif strfind(PARM{n,1},'COMPASS_')
        elseif strfind(PARM{n,1},'RC')
        elseif strfind(PARM{n,1},'SERVO')
        elseif strfind(PARM{n,1},'BATT_')
        elseif strfind(PARM{n,1},'AHRS_')
        elseif strfind(PARM{n,1},'FS_')
        %elseif strfind(PARM{n,1},'Q_')

        elseif isnan(index)
            fprintf('%s, Does not exist!!\n',PARM{n,1});
        elseif PARM{n,2} ~= PRAM2.PARM{index,2}
            fprintf('%s, %g, %g\n',PARM{n,1},PARM{n,2},PRAM2.PARM{index,2});
        end
        
    end
    return
end






% Caculaate P
%{
% Parameters
P_gain = 2;
I_gain = 0.6;
D_gain = 0.1;
scaleing_speed = 15;
trim_throttle = 45;
R_max = 60;
rll_t_const = 0.5;

P_gain = (P_gain - I_gain*rll_t_const) * rll_t_const - D_gain;

index1 = find((ATT(:,2)/1e6)<TimeCrop(1),1,'last');
index2 = find((ATT(:,2)/1e6)>TimeCrop(2),1,'first');

for n = index1:index2
  
    % interpolate airspeed time to attitude time
    airspd = interp1(CTUN(:,2),CTUN(:,10),ATT(n,2));
    throttle = interp1(AETR(:,2),AETR(:,5),ATT(n,2));
       
    Roll_error = ATT(n,3)-ATT(n,4);
    
    Roll_error = Roll_error * (1/rll_t_const);
    
    if Roll_error > R_max
        Roll_error = R_max;
    elseif Roll_error < -R_max
        Roll_error = -R_max;
    end

    
    if isnan(airspd)
        fprintf('NaN Aspd @ t = %g s\n',ATT(n,2)/1e6)
    end
    
    % Throttle scaling
%{
    if throttle > 0
        scalling = 0.5 + (trim_throttle / throttle / 2); % First order taylor expansion of square root
    else
        scalling = 1.67;
    end
    if scalling > 1.67
        scalling = 1.67;
    elseif scalling < 0.6
        scalling = 0.6;
    end
%}
    
    % Airspeed scailing
    scalling = scaleing_speed/airspd;
    
    if scalling > 2
        scalling = 2;
    elseif scalling < 0.5
        scalling = 0.5;
    end
    
    P_calc((n-index1)+1) = Roll_error * P_gain * scalling;
    P_calc2((n-index1)+1) = -R_max * P_gain * scalling;
    sailing_calc((n-index1)+1) = scalling;
end
plot(ATT((index1:index2),2)/1e6,P_calc)
plot(ATT((index1:index2),2)/1e6,P_calc2)

figure
hold all
plot(ATT((index1:index2),2)/1e6,sailing_calc)
if ~isempty(TimeMark) && Time_markers == 1
    for n = 1:length(TimeMark)
        %Plotting Time Event Markers
        plot([TimeMark(n),TimeMark(n)],[0.5,2],':g')
    end
end
%ylim([-45,45])
return
%}
%}

% try and load in defualt file
if exist(file_name,'file')
    Load_file([],[],file_name)
end

end

% select a log file
function select_file(~,~)
[file,path] = uigetfile('*.mat','Pick Log File');

handles = getappdata(0,'handles');

% only update file name if has been selected
if file ~= 0
    handles.open_param_file_default.String = file;
    
    setappdata(0,'file',[path,file]);
    
    Load_file();
end

end

% load in a log file
function Load_file(~,~,file_name)
if isappdata(0,'file')
    file_name_user = getappdata(0,'file');
    if exist(file_name_user,'file')
        file_name = file_name_user;
    end
end

if ~exist(file_name,'file')
    fprintf('Could not find file: %s\n\n',file_name)
    return
end

clc
logs = load_log(file_name);
setappdata(0,'logs',logs);

handles = getappdata(0,'handles');

% Print messages
MSG1 = get_var('MSG1');
if strcmp(MSG1,'false') && handles.print_messages.Value == 1
    fprintf('Could not find messages\n\n')
else
    for n = 1:length(MSG1.value)
        message = MSG1.value{1,n}{3};
        if handles.print_messages.Value == 1
            fprintf('%gs - %s\n',MSG1.value{1,n}{2}/1e6, message)
        end
        % Try and spot if copter plane or rover and set window name
        if contains(message,'Copter')
            setappdata(0,'vehicle','Copter');
            set(handles.figure,'Name',message,'NumberTitle','off')
        elseif contains(message,'Plane')
            setappdata(0,'vehicle','Plane');
            set(handles.figure,'Name',message,'NumberTitle','off')            
        elseif contains(message,'Rover')
            setappdata(0,'vehicle','Rover');
            set(handles.figure,'Name',message,'NumberTitle','off')            
        end
    end
    if handles.print_messages.Value == 1
        fprintf('\n')
    end
end

% save to  file
txt_file_name = [erase(file_name,'.mat'),'.param'];
if handles.save_param_file.Value == 1
    PARM = get_var('PARM');
    if strcmp(PARM,'false')
        fprintf('Could not find params\n\n')
    else
        fid = fopen(txt_file_name, 'w+');
        for n = 1:length(PARM.value)
            fprintf(fid,'%16s\t%g\n',PARM.value{n,1},PARM.value{n,2});
        end
        fclose(fid);
    end
end

% open file
if handles.open_param_file.Value == 1
    if exist(txt_file_name,'file')
        open(txt_file_name)
    else
        fprintf('Error, could not find pramiter file: %s\n\n',txt_file_name)
    end
end

% ungrey out relevent buttons
handles.mode_shadeing.Enable = 'on';
handles.time_markers.Enable = 'on';
handles.time_crop.Enable = 'on';

plots = getappdata(0,'plots');
plot_catagorys = getappdata(0,'plot_catagorys');

% check for plot button avalablity
% avalable = zeros(length(plots),1);
plot_names = cell(length(plots),1);
for n = 1:length(plots)
    i = 1;
    while i <= (length(plots{n})-4)/2
        avalable(n,i) = ~strcmp('false', get_var(plots{n}{3+(i*2)}));%#ok<AGROW>
        i = i + 1;
    end
    plot_names{n} = plots{n}{1};
end

% any plot from each catogorcy can be avalable
for n = 1:length(plot_catagorys)
    index = find(strcmp(plot_names,plot_catagorys{n}));
    if any(any(avalable(index,:)))
        handles.(plots{index(1)}{1}).Enable = 'on';
    end
end

if any(any(avalable))
    handles.plot.Enable = 'on';
end

% Save avalibility array to save redoing
setappdata(0,'avalable',avalable);

end

function plot_figures(~,~)
plot_close();

plots = getappdata(0,'plots');
%plot_catagorys = getappdata(0,'plot_catagorys');
handles = getappdata(0,'handles');
avalable = getappdata(0,'avalable');

% plot_log('ESC Log','TEMP (deg C)',{ESC1,ESC1_label},6,{ESC2,ESC2_label},6)


% Plot all the selected plots 
for n = 1:length(plots)
    if any(avalable(n,:)) && handles.(plots{n}{1}).Value
        plot_vars = [];
        plot_vars{1} = plots{n}{3};
        plot_vars{2} = plots{n}{4};
        plot_vars{3} = plots{n}{2};
        for i = 5:2:3+(sum(avalable(n,:))*2)
           plot_vars{i-1} =  get_var({plots{n}{i}});
           plot_vars{i} = plots{n}{i+1};
        end
        
        for i = length(plot_vars):-1:4
            if strcmp(plot_vars{i},'false')
                plot_vars(i+1) = [];
                plot_vars(i) = [];
            end
        end
        
        if length(plot_vars) ~= 3
            plot_log(plot_vars{:})
        else
            fprintf('Could not find Data for plot: %s\n',plot_vars{1})
        end
    end
end

handles.plot_close.Enable = 'on';
end

function plot_close(~,~)
% Close figs all except GUI window
all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, 1));

handles = getappdata(0,'handles');
handles.plot_close.Enable = 'off';
end

% get data for a particular varable
function variable = get_var(var_name)
logs = getappdata(0,'logs');

for n = 1:length(logs)
    if strcmp(logs{n,1},var_name)
        index = n;
        break;
    end
end

if exist('index','var')
    variable.name = logs{index,1};
    variable.label = logs{index,2};
    variable.value = logs{index,3};
else
    variable = 'false';
end
end

% Plot function
function plot_log(Title,Label,group,varargin)

colour = ['b','r','m','g','k','c','b','r','m','g','k','c','b','r','m','g','k','c','b','r','m','g','k','c','b','r','m','g','k','c','b','r','m','g','k','c','b','r','m','g','k','c','b','r','m','g','k','c','b','r','m','g','k','c'];
% Mdes for plane
% Mode for copter

vehicle = getappdata(0,'vehicle');
if strcmp(vehicle,'Copter')
    MODELOOKUP = {'Stabilize','Acro','AltHold','Auto','Guided','Loiter','RTL','Circle','Land','Drift','Sport','Flip','AutoTune','PosHold','Brake','Throw','Avoid_ADSB','Guided_NoGPS','Smart_RTL','FlowHold','Follow'};
elseif strcmp(vehicle,'Plane')
    MODELOOKUP = {'Manual','CIRCLE','STABILIZE','TRAINING','ACRO','FBWA','FBWB','CRUISE','AUTOTUNE','','Auto','RTL','Loiter','','AVOID_ADSB','Guided','','QSTABILIZE','QHOVER','QLOITER','QLAND','QRTL'};
elseif  strcmp(vehicle,'Rover')
    MODELOOKUP = {'Not Done this yet'};
end

figure('name',group,'NumberTitle','off')
hold all

legend_val = {};

Min_max_x = [inf,-inf];
Min_max_y = [inf,-inf];

for n = 1:(length(varargin)/2)
    var{1} = varargin{(n*2)-1}.value;
    var{2} = varargin{(n*2)-1}.label;

    var_index = varargin{(n*2)};
    
    for m = 1:length(var_index)
        if sign(var_index(m)) == 1
            
            plot(var{1}(:,2)/1e6,var{1}(:,var_index(m)),colour(var_index(m) + n - 1))
            
            legend_val(end+1) = var{2}(var_index(m));
            
            Min_max_x = [min([Min_max_x(1),min(var{1}(:,2)/1e6)]),max([Min_max_x(2),max(var{1}(:,2)/1e6)])];
            Min_max_y = [min([Min_max_y(1),min(var{1}(:,var_index(m)))]),max([Min_max_y(2),max(var{1}(:,var_index(m)))])];
        else
            var_index(m) = abs(var_index(m));
            
            plot(var{1}(:,2)/1e6,-1 * var{1}(:,var_index(m)),colour(m))
            
            legend_val(end+1) = strcat('NEGATIVE ',var{2}(var_index(m)));
            
            Min_max_x = [min([Min_max_x(1),min(var{1}(:,2)/1e6)]),max([Min_max_x(2),max(var{1}(:,2)/1e6)])];
            Min_max_y = [min([Min_max_y(1),-1*max(var{1}(:,var_index(m)))]),max([Min_max_y(2),-1*min(var{1}(:,var_index(m)))])];
        end
        
        if any(isnan(var{1}(:,var_index(m))))
            if min(isnan(var{1}(:,var_index(m))))
                legend_val(end) = strcat(legend_val(end),'\color{red} Completely NAN!');
            else
                legend_val(end) = strcat(legend_val(end),'\color{red} Contains NAN!');
            end
        end
        
    end
    
end

%{
%--- Manual Time Markers ---
if ~isempty(TimeMark) && Time_markers == 1
    for n = 1:length(TimeMark)
        %Plotting Time Event Markers
        plot([TimeMark(n),TimeMark(n)],Min_max_y,':g')
    end
end

%--- Mode Switch Markers ---
if Mode_markes
    for n = 1:size(MODE,1)
        %Plotting Time Event Markers
        plot([MODE(n,2)/1e6,MODE(n,2)/1e6],Min_max_y,'-.k')
        text((MODE(n,2)/1e6),Min_max_y(1)*0.9,MODELOOKUP{MODE(n,4)+1})
    end
end
%}
xlabel('Time (s)')

ylabel(Label)
title(Title)
legend(legend_val)
%{
if sum(isnan(TimeCrop)) == 2
    xlim(Min_max_x)
    ylim(Min_max_y)
elseif isnan(TimeCrop(1))
    xlim([Min_max_x(1),TimeCrop(2)])
elseif isnan(TimeCrop(2))
    xlim([TimeCrop(1),Min_max_x(2)])
else
    xlim(TimeCrop)
end
%}
end

% load in a log
function logs = load_log(file_name)

load(file_name);

vars = who;
i = 1;
for n = 1:length(vars)
    name = vars{n};
    
    % All varables posible must have _label
    if endsWith(name, '_label')
        % some varables might have data aswell
        name = erase(name,'_label');
        if any(strcmp(vars,name))
            logs{i,1} = name;
            logs{i,2} = eval([name,'_label']);
            logs{i,3} = eval(name);
            i = i + 1;
        end
    end
    
    % except MSG1 and PARM, and probibly some others I havent noticed
    if strcmp(name,'MSG1') || strcmp(name,'PARM')
        logs{i,1} = name;
        logs{i,3} = eval(name);
        i = i + 1;
    end
end
end
