function Log_Plotter ()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
clc
clear
close all
if isappdata(0,'file')
    rmappdata(0,'file');
end
if isappdata(0,'time_crop')
    rmappdata(0,'time_crop');
end
if isappdata(0,'time_markers')
    rmappdata(0,'time_markers');
end


% posible plots to pick from
% button varable name, button display name, Plot title, Plot Y axis label, Plot varable, Plot colums in varable, varable 2, colums in varable 2 ect
plots{1} = {'att','Attitude','Roll Angle - Desired vs actual','Angle (deg)','ATT',[3,4]};
plots{2} = {'att','Attitude','Pitch Angle - Desired vs actual','Angle (deg)','ATT',[5,6]};
plots{3} = {'att','Attitude','Yaw Angle - Desired vs actual','Angle (deg)','ATT',[7,8]};
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
plots{14}= {'Speed','Speed','Speed (m/s)','GPS',11,'GPS2',11,'CTUN',10};
plots{15}= {'PWM_Out','PWM Out','PWM Out','PWM','RCOU',[3,4,5,6,7,8,9,10]};
plots{16}= {'Control_outputs','Control out','Control outputs','','AETR',[3,4,5,6,7]};
plots{17}= {'Quadplane_PID','Quadplane PID','Quadplane Pitch PID','','PIQP',[4,5,6]};
plots{18}= {'Quadplane_PID','Quadplane PID','Quadplane Roll PID','','PIQR',[4,5,6]};
plots{19}= {'Quadplane_PID','Quadplane PID','Quadplane Yaw PID','','PIQY',[4,5,6]};
plots{20}= {'AoA','Angle of Attack','Angle of Attack','AoA (deg)','AOA',3};
plots{21}= {'Nav','Navigation','Navigation','','NTUN',6};
plots{22}= {'Bat','Battery','Battery Volts','Voltage','BAT',4};
plots{23}= {'Bat','Battery','Battery Current','Current (A)','BAT',5};
plots{24}= {'rate','Rates','Roll Rate','rate deg/s^2','RATE',[3,4]};
plots{25}= {'rate','Rates','Pitch Rate','rate deg/s^2','RATE',[6,7]};
plots{26}= {'rate','Rates','Yaw Rate','rate deg/s^2','RATE',[9,10]};
plots{27}= {'ESC_Log','ESC Log','ESC RPM','RPM','ESC1',3,'ESC2',3,'ESC3',3,'ESC4',3,'ESC5',3,'ESC6',3,'ESC7',3,'ESC8',3};
plots{28}= {'ESC_Log','ESC Log','ESC Voltage','Voltage','ESC1',4,'ESC2',4,'ESC3',4,'ESC4',4,'ESC5',4,'ESC6',4,'ESC7',4,'ESC8',4};
plots{29}= {'ESC_Log','ESC Log','ESC Current','Current (A)','ESC1',5,'ESC2',5,'ESC3',5,'ESC4',5,'ESC5',5,'ESC6',5,'ESC7',5,'ESC8',5};
plots{30}= {'ESC_Log','ESC Log','ESC Tmep','Temp (deg C)','ESC1',6,'ESC2',6,'ESC3',6,'ESC4',6,'ESC5',6,'ESC6',6,'ESC7',6,'ESC8',6};
plots{31}= {'RC_contorl_in','Control in2','Control inputst','RCIN','RCIN',[3,4,5,6,7,8,9,10,11,12,13,14,15,16]};


% Defualts
file_name = 'Input File Name';
print_messages_default = 1;
save_param_file_default = 1;
open_param_file_default = 0;
mode_shading_default = 1;
time_markers_default = 0;
time_crop_default = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs if you dont want to use GUI
% copy varable from defaults and over write


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple(ish) GIU
load_in = @(src,event)Load_file(src,event,file_name);

handles.figure = figure('MenuBar','none','Resize', 'off');


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
        'Position',[X(n)+10 (Y(n)+y_spacing*(Y_no+1)+15) 100 30],...
        'Enable','off',...
        'Value',0,...
        'HandleVisibility','off');
end

fig_hight = y_spacing*(Y_no+1)+15; 

handles.print_messages = uicontrol('Style',...
    'checkbox',...
    'String','Print MSG''s',...
    'Position',[10 fig_hight+(y_spacing*3)+10 100 30],...
    'Value',print_messages_default,...
    'TooltipString','Print messages to comand window',...
    'HandleVisibility','off');

handles.save_param_file = uicontrol('Style',...
    'checkbox',...
    'String','Save Param File',...
    'Position',[110 fig_hight+(y_spacing*3)+10 100 30],...
    'Value',save_param_file_default,...
    'TooltipString','Save parmiters from log to .param file',...
    'HandleVisibility','off');

handles.open_param_file = uicontrol('Style',...
    'checkbox',...
    'String','Open Param File',...
    'Position',[210 fig_hight+(y_spacing*3)+10 100 30],...
    'Value',open_param_file_default,...
    'TooltipString',['Open .param with the same name as current log' newline 'file can be created using save param file'],...
    'HandleVisibility','off');

handles.open_param_file_default = uicontrol('Style',...
    'pushbutton',...
    'String',file_name,...
    'Position',[10 fig_hight+(y_spacing*2)+5 200 30],...
    'HandleVisibility','off',...
    'Callback', @select_file);

handles.load_file = uicontrol('Style',...
    'pushbutton',...
    'String','Reload',...
    'Position',[210 fig_hight+(y_spacing*2)+5 100 30],...
    'HandleVisibility','off',...
    'Callback', load_in);

handles.mode_shadeing = uicontrol('Style',...
    'checkbox',...
    'String','mode shading',...
    'Position',[10 fig_hight+y_spacing 100 30],...
    'Value',mode_shading_default,...
    'Enable','off',...
    'TooltipString','Shade plot backgrounds with coulors representing modes and generate key',...
    'HandleVisibility','off');

handles.time_markers = uicontrol('Style',...
    'checkbox',...
    'String','time markers',...
    'Position',[110 fig_hight+y_spacing 100 30],...
    'Value',time_markers_default,...
    'Enable','off',...
    'TooltipString','Plot vertical lines at specified times',...
    'HandleVisibility','off',...
    'Callback', @get_time_markers);

handles.time_crop = uicontrol('Style',...
    'checkbox',...
    'String','time crop',...
    'Position',[210 fig_hight+y_spacing 100 30],...
    'Value',time_crop_default,...
    'Enable','off',...
    'TooltipString','Crop plot to specified time period',...
    'HandleVisibility','off',...
    'Callback', @get_time_crop);

handles.plot = uicontrol('Style',...
    'pushbutton',...
    'String','Plot',...
    'Position',[110 10 100 30],...
    'Enable','off',...
    'HandleVisibility','off',...
    'Callback', @plot_figures);

handles.plot_close = uicontrol('Style',...
    'pushbutton',...
    'String','Close Plots',...
    'Position',[10 10 100 30],...
    'Enable','off',...
    'HandleVisibility','off',...
    'Callback', @plot_close);

position = handles.figure.Position;
position(3) = 320;  
position(4) = fig_hight+(y_spacing*4)+15;
handles.figure.Position = position;

setappdata(0,'handles', handles);
setappdata(0,'plots', plots);
setappdata(0,'plot_catagorys', plot_catagorys);

% try and load in defualt file
if exist(file_name,'file')
    Load_file([],[],file_name)
end

end

function get_time_crop(~,~)
handles = getappdata(0,'handles');

if handles.time_crop.Value == 0
    return;
end

if isappdata(0,'time_crop')
    default_time_crop = getappdata(0,'time_crop');
else
    default_time_crop = {'',''};
end

time_crop = inputdlg({'Start time (s)','End time (s)'},'Time Crop',[1 35],default_time_crop);

if isempty(time_crop(:))
    return;
end
    test = str2num(time_crop{1}); %#ok<ST2NM>
    if isempty(test)
        time_crop{1} = '';
    else
        time_crop{1} = num2str(test(1));
    end
    
    test = str2num(time_crop{2}); %#ok<ST2NM>
    if isempty(test)
        time_crop{2} = '';
    else
        time_crop{2} = num2str(test(1));
    end
    
    if str2num(time_crop{2}) <= str2num(time_crop{1}) %#ok<ST2NM>
        fprintf('Time crop - Cannot have start after end\n')
        % if only changed one then reset it 
        if str2num(time_crop{2}) == str2num(default_time_crop{2}) %#ok<ST2NM>
            time_crop{1} = '';
        else
            time_crop{2} = '';           
        end
    end
    
    setappdata(0,'time_crop',time_crop);    
end

function get_time_markers(~,~)
handles = getappdata(0,'handles');

if handles.time_markers.Value == 0
    return;
end

if isappdata(0,'time_markers')
    default_time_mark = getappdata(0,'time_markers');
else
    default_time_mark = {''};
end

time_mark = inputdlg('Time Markers(s) seperate with spaces','Time Markers',[1 40],default_time_mark);

if isempty(time_mark(:))
    return;
end
    
    test = str2num(time_mark{1}); %#ok<ST2NM>
    if isempty(test)
        time_mark{1} = '';
    end    
    setappdata(0,'time_markers',time_mark);
    
end

% select a log file
function select_file(~,~)
if isappdata(0,'file')
    file = getappdata(0,'file');
    [file,path] = uigetfile('*.mat','Pick Log File',file);
else
    [file,path] = uigetfile('*.mat','Pick Log File');
end


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

% Print any error's 
ERR = get_var('ERR');
errors = {'MAIN','RADIO','COMPASS','OPTFLOW','FAILSAFE RADIO','FAILSAFE BATT','FAILSAFE GPS','FAILSAFE GCS','FAILSAFE FENCE','FLIGHT MODE','GPS','CRASH CHECK','FLIP','AUTOTUNE','PARACHUTE','EKFCHECK','FAILSAFE EKFINAV','BARO','CPU','FAILSAFE ADSB','TERRAIN','NAVIGATION','FAILSAFE TERRAIN','EKF PRIMARY','THRUST LOSS CHECK'};
sub_code = [4,    2,      NaN,      NaN,       3,               3,              3,             3,             3,              NaN,           12,   5,            6,    NaN       9,          10,        NaN,                11,   NaN,  3,              7,        8,           3,                 NaN,           NaN];

% general error codes
errors_sub{1+0,1} = 'ERROR RESOLVED';
errors_sub{1+1,1} = 'FAILED TO INITIALISE';
errors_sub{1+4,1} = 'UNHEALTHY';

% radio
errors_sub{1+2,2} = 'LATE FRAME';

% failsafe_thr, batt, gps
errors_sub{1+0,3} = 'FAILSAFE RESOLVED';
errors_sub{1+1,3} = 'FAILSAFE OCCURRED';

% main
errors_sub{1+1,4}  = 'INS DELAY';

% crash checker
errors_sub{1+1,5}  = 'CRASH';
errors_sub{1+2,5}  = 'LOSS OF CONTROL';

% flip
errors_sub{1+2,6}  = 'ABANDONED';

% terrain
errors_sub{1+2,7} = 'MISSING TERRAIN DATA';

% navigation
errors_sub{1+2,8} = 'FAILED TO SET DESTINATION';
errors_sub{1+3,8} = 'RESTARTED RTL';
errors_sub{1+4,8} = 'FAILED CIRCLE INIT';
errors_sub{1+5,8} = 'DEST OUTSIDE FENCE';

% parachute
errors_sub{1+2,9} = 'TOO LOW';
errors_sub{1+3,9} = 'LANDED';

% EKF check definitions
errors_sub{1+2,10} = 'BAD VARIANCE';
errors_sub{1+0,10} = 'VARIANCE CLEARED';

% Baro specific error codes
errors_sub{1+2,11} = 'GLITCH';

% GPS specific error coces
errors_sub{1+2,12} = 'GLITCH';

if ~strcmp(ERR,'false')
    error_time = zeros(1,length(ERR.value));
    error_print = cell(1,length(ERR.value));
    for n = 1:length(ERR.value)
        error_time(n) = ERR.value(n,2)/1e6;
        message = errors(ERR.value(n,3));
        index = sub_code(ERR.value(n,3));
        index2 = ERR.value(n,4) + 1;
        if isnan(index) || isempty(errors_sub{index,index2})
            index = 1;
        end
        mesage_2 = errors_sub(index,index2);
        error_print{n} = sprintf('%gs - %s %s\n',error_time(n), message{:},mesage_2{:});
    end
end

% Print messages
MSG1 = get_var('MSG1');
if strcmp(MSG1,'false') && handles.print_messages.Value == 1
    fprintf('Could not find messages\n\n')
    if exist('error_time','var')
        for n = 1:length(error_time)
            if  exist('cprintf','file') == 2
                cprintf('red','%s',error_print{n})
            else
                fprintf('%s',error_print{n})
            end
        end
    end
else
    i = 1;
    for n = 1:length(MSG1.value)
        message = MSG1.value{1,n}{3};
        if handles.print_messages.Value == 1
            time = MSG1.value{1,n}{2}/1e6;
            
            if exist('error_time','var')
                while time > error_time(i)
                    if  exist('cprintf','file') == 2
                        cprintf('red','%s',error_print{i})
                    else
                        fprintf('%s',error_print{i})
                    end
                    i = i + 1;
                end
            end
            fprintf('%gs - %s\n',time, message)
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
        elseif contains(message,'AntennaTracker')
            setappdata(0,'vehicle','AntennaTracker');
            set(handles.figure,'Name',message,'NumberTitle','off')
        end
    end
    
end
if exist('error_time','var')
    for n = i:length(error_time)
        if  exist('cprintf','file') == 2
            cprintf('red','%s',error_print{n})
        else
            fprintf('%s',error_print{n})
        end
    end
end
if handles.print_messages.Value == 1
    fprintf('\n')
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

vehicle = getappdata(0,'vehicle');
handles = getappdata(0,'handles');

if strcmp(vehicle,'Copter')
    MODELOOKUP = {'Stabilize','Acro','AltHold','Auto','Guided','Loiter','RTL','Circle','Land','Drift','Sport','Flip','AutoTune','PosHold','Brake','Throw','Avoid_ADSB','Guided_NoGPS','Smart_RTL','FlowHold','Follow'};
elseif strcmp(vehicle,'Plane')
    MODELOOKUP = {'Manual','CIRCLE','STABILIZE','TRAINING','ACRO','FBWA','FBWB','CRUISE','AUTOTUNE','','Auto','RTL','Loiter','','AVOID_ADSB','Guided','','QSTABILIZE','QHOVER','QLOITER','QLAND','QRTL'};
elseif  strcmp(vehicle,'Rover')
    MODELOOKUP = {'Not Done this yet'};
elseif  strcmp(vehicle,'AntennaTracker')
    MODELOOKUP = {'MANUAL','STOP','SCAN','SERVO TEST'};
    MODELOOKUP{11} = 'AUTO';
    MODELOOKUP{17} = 'INITIALISING';
end

fig_handle = figure('name',group,'NumberTitle','off');
position = fig_handle.Position;
position(3) = position(3) * 1.25;  
fig_handle.Position = position;

ax = axes;
hold all

legend_val = {};

Min_max_x = [inf,-inf];
Min_max_y = [inf,-inf];

for n = 1:(length(varargin)/2)
    var{1} = varargin{(n*2)-1}.value;
    var{2} = varargin{(n*2)-1}.label;
    var{3} = varargin{(n*2)-1}.name;

    var_index = varargin{(n*2)};
    
    legend_val = cell(length(var_index),1);
    for m = 1:length(var_index)
        if sign(var_index(m)) == 1
            
            plot(ax,var{1}(:,2)/1e6,var{1}(:,var_index(m)),colour(var_index(m) + n - 1))
            
            legend_val(m) = strcat(var{3},{' - '},var{2}(var_index(m)));
            
            Min_max_x = [min([Min_max_x(1),min(var{1}(:,2)/1e6)]),max([Min_max_x(2),max(var{1}(:,2)/1e6)])];
            Min_max_y = [min([Min_max_y(1),min(var{1}(:,var_index(m)))]),max([Min_max_y(2),max(var{1}(:,var_index(m)))])];
        else
            var_index(m) = abs(var_index(m));
            
            plot(ax,var{1}(:,2)/1e6,-1 * var{1}(:,var_index(m)),colour(m))
            
            legend_val(m) = strcat(var{3},{' - '},{'NEGATIVE '},var{2}(var_index(m)));
            
            Min_max_x = [min([Min_max_x(1),min(var{1}(:,2)/1e6)]),max([Min_max_x(2),max(var{1}(:,2)/1e6)])];
            Min_max_y = [min([Min_max_y(1),-1*max(var{1}(:,var_index(m)))]),max([Min_max_y(2),-1*min(var{1}(:,var_index(m)))])];
        end
        
        if any(isnan(var{1}(:,var_index(m))))
            if min(isnan(var{1}(:,var_index(m))))
                legend_val(m) = strcat(legend_val(m),'\color{red} Completely NAN!');
            else
                legend_val(m) = strcat(legend_val(m),'\color{red} Contains NAN!');
            end
        elseif all((var{1}(:,var_index(m))) == 0)
            legend_val(m) = strcat(legend_val(m),'\color{red} all zero!');
        end
        
        
    end
    
end
if any(isinf(Min_max_y)) || all(Min_max_y == 0)
    Min_max_y = [0,1];
end


%--- Manual Time Markers ---BCL4_label
if handles.time_markers.Value == 1 && isappdata(0,'time_markers')    
    stringin = getappdata(0,'time_markers');
    TimeMark = str2num(stringin{1}); %#ok<ST2NM>
    for n = 1:length(TimeMark)
        %Plotting Time Event Markers
        plot(ax,[TimeMark(n),TimeMark(n)],Min_max_y,':g','linewidth',2)
    end
end

%--- Time crop ---
if handles.time_crop.Value == 1 && isappdata(0,'time_crop')
    stringin = getappdata(0,'time_crop');
    test = str2num(stringin{1}); %#ok<ST2NM>
    if isempty(test)
        TimeCrop(1) = Min_max_x(1);
    else
        TimeCrop(1) = test;
    end
    
    test = str2num(stringin{2}); %#ok<ST2NM>
    if isempty(test)
        TimeCrop(2) = Min_max_x(2);
    else
        TimeCrop(2) = test;
    end
    Min_max_x = TimeCrop;
end

xlim(ax,Min_max_x)
ylim(ax,Min_max_y)
xlabel(ax,'Time (s)')
ylabel(ax,Label)
title(ax,Title)
legend(ax,legend_val,'Location','northeastoutside');

%--- Mode Switch Markers ---
if handles.mode_shadeing.Value == 1
    varable = get_var('MODE');
    if ~strcmp(varable,'false')
        ax2 = axes('visible','off');
        
        MODE = varable.value;
        
        %Appending final timestep onto end of mode marker lookup
        ModeTime = [MODE(:,2)/1e6;Min_max_x(2)];
        
        C = colormap('hsv');
        
        modes = unique(MODE(:,4)+1);
        patch_hand = zeros(size(MODE,1) ,1);
        legend_val_2 = cell(size(MODE,1) ,1);
        for n = 1:size(MODE,1) 
            
            colour_val = round((find(MODE(n,4)+1 == modes) / length(modes)) * 64);

            %Plotting time event markers as colour patches
            patch_hand(n) = patch([ModeTime(n),ModeTime(n),ModeTime(n+1),ModeTime(n+1),ModeTime(n)],...
                [Min_max_y(1),Min_max_y(2),Min_max_y(2),Min_max_y(1),Min_max_y(1)],...
                C(colour_val,:),...
                'FaceAlpha',0.15,...
                'EdgeColor','none');
                       
            legend_val_2{n} = MODELOOKUP{1,MODE(n,4)+1};
            
        end
        [~, index] = unique(legend_val_2,'stable');
        legend(ax2,patch_hand(index),legend_val_2{index},'Location','southeastoutside')
        xlim(ax2,Min_max_x)
        ylim(ax2,Min_max_y)
        set(ax2,'position',get(ax,'position'))
    end
end
end

% load in a log
function logs = load_log(file_name)

load(file_name); %#ok<LOAD>

vars = who;
i = 1;
for n = 1:length(vars)
    name = vars{n};
    
    % All varables posible must have _label
    if endsWith(name, '_label')
        % some varables might have data aswell
        name = erase(name,'_label');
        if any(strcmp(vars,name))
            logs{i,1} = name; %#ok<AGROW>
            logs{i,2} = eval([name,'_label']); %#ok<AGROW>
            logs{i,3} = eval(name); %#ok<AGROW>
            i = i + 1;
        end
    end
    
    % except MSG1 and PARM, and probibly some others I havent noticed
    if strcmp(name,'MSG1') || strcmp(name,'PARM')
        logs{i,1} = name; %#ok<AGROW>
        logs{i,3} = eval(name); %#ok<AGROW>
        i = i + 1;
    end
end
end
