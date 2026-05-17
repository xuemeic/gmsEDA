%%%%%%% get data
load("eda_5m32s_23ptcp_shift.mat")
% this loads a structure with field eda and id
eda = eda_5m32s_23ptcp_shift.eda;
id_all = eda_5m32s_23ptcp_shift.id;

eda_len = size(eda, 1);

%%%%%% process one signal at a time
idx = 2; % pick 1 through 23
signal = eda(:,idx);

[phasic, baseline] = gmsEDA(signal);

%%%%%%%%% plot
uncw_blue = [15, 44, 88]/255;
total_sec = 5*60 + 32;
timeshift = 6; % in seconds. 

xt = 0:30:(total_sec-timeshift); % xticks
x_time = linspace(0, eda_len/4, eda_len); % in seconds

f = figure(1);
f.Position = [600, 600, 700, 500];
subplot(2,1,1)
plot(x_time(1:(eda_len-timeshift*4)), eda((timeshift*4+1):eda_len, idx), 'Color',uncw_blue)
id = strsplit(id_all{idx}, '_'); % 'SUD_003_14_EDA'
% the third one gives the number
id = id{3};
t = title(['Subject ', id]);
t.FontSize = 15;
yt = ylabel('raw EDA');
yt.FontSize = 15;
clear xticks
xlim([0, total_sec-timeshift])
xticks(xt);
hold on

subplot(2,1,2)
plot(x_time(1:(eda_len-timeshift*4)), phasic((timeshift*4+1):eda_len), 'Color',uncw_blue);
yt2 = ylabel('phasic');
yt2.FontSize = 15;
xlim([0, total_sec-timeshift])
xticks(xt);

