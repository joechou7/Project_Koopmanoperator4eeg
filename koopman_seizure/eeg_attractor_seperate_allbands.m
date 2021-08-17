addpath('./utils');
%

path = '../DATA/chb-mit-scalp-eeg/';
d = dir(path);

dt = 1/256;
windows_second = 50;
time_window = windows_second*256;
stackmax = 500;  % number of shift-stacked rows
rmax = 15;
overlap = 0.8;
seizure_percent = 0.7;
names = ["filtdelta","filttheta","filtalpha","filtbeta","filtgamma"];
sig_mean = 0;
method = 'max';
numoffeature = 23*6*5;

load('channels_map.mat')

for i = 12+7 %13:(length(d)-1)
    foldername = d(i).name;
    %     di = dir([path,foldername]);
    fid = fopen([path,foldername,'/',foldername,'-summary.txt']);
    tline = fgetl(fid);
    sessions = [];
    index = [];
    while ischar(tline)
        if contains(tline, 'Channels ')
            tline = fgetl(fid);
            tline = fgetl(fid);
            blank = 0;
            index = [];
            first_d = true;
            while contains(tline, 'Channel ')
                num = regexp(tline,'Channel \d+:','match');
                num = num{1};
                num = str2double(num(9:end-1))-blank;
                name = regexp(tline,': .*','match');
                name = name{1}(3:end);
                if strcmp(name,'T8-P8')
                    if ~first_d 
                        name = [name,'-2'];
                    end
                    first_d = false;
                end
                if isKey(channels, name)
                    num_c = channels(name);
                    temp = [num;num_c];
                    index = [index,temp];
                elseif strcmp(name,'-')
                    blank = blank+1;
                end
                
                tline = fgetl(fid);
            end
            tline = fgetl(fid);    
        end
        if size(index,2) ~= 23
            tline = fgetl(fid);
            continue
        end
        
        if contains(tline, 'File Name')
            session = regexp(tline,'_\d+.','match');
            filename = [tline(12:end-4),'.mat'];
            session = session{1};
            session = str2double(session(2:end-1));
            if session < 1
                tline = fgetl(fid);
                continue
            end
        else
            tline = fgetl(fid);
            continue
        end
        while ~contains(tline, 'Number of Seizures in File')
            tline = fgetl(fid);
        end
        num = str2double(tline(29));
        if num == 0
            continue
        end
        
        % load data
        load([path, foldername,'/',filename]);
        data = data.eeg;
        
        data = data(~isnan(data(:,1)),:);
        
        [~, order] = sort(index(2,:));
        index = index(:,order);
        
        data = data(index(1,:)',:);
        
        [winstart,nw,mod_end,interval] = slicewin(size(data,2),time_window,overlap);
        sessions = nw;
        % read file obtain seizure timestamp
        label = zeros(nw-1,1);
        index_seizure = [];
        start = -1;
        for nextline = 1:num
            tline = fgetl(fid);
            start = regexp(tline,'\d+ second','match');
            start = start{1};
            start = str2double(start(1:end-7));
            tline = fgetl(fid);
            ends = regexp(tline,' \d+ second','match');
            ends = ends{1};
            ends = str2double(ends(2:end-7));
            index_seizure = label_seizure(start/dt,ends/dt,interval,seizure_percent);
            label(index_seizure) = 1;
        end
        
        
        xdata = zeros(size(data,1),time_window,nw-1);
        for n = 1:nw
            xdata(:,:,n) = data(:,winstart(n):(winstart(n)+time_window-1));
        end
        
        res = zeros(nw-1,numoffeature);
        for window = 1:nw-1
            temp = [];
            for channel = 1:size(xdata,1)
                data_t = xdata(channel,:,window);
                %                 data_t = data(channel,:);
                [cleandata,filtdelta,filttheta,filtalpha,filtbeta,filtgamma] = eeg_bands(data_t,3);
                
                
                
                temp_data = filtdelta;
                [U,S,V,r,tspan,x,y,t] = hankle_analysis_v2(temp_data',dt, stackmax, rmax);
                %                 figure;plot3(V(:,1),V(:,2),V(:,3))
                V = V(:,1:r);
                %                 sig_mean = mean(V,1);
                ev = diag(S);
                ev = ev(1:r,1);
                ev = ev./sum(ev);
                %                 tt_delta = sqrt(((V-sig_mean).^2)*ev);
                tt_delta = system_recon(sig_mean,V,ev,method);
                %                 figure;plot(tt_delta)
                %                 figure;plot(system_recon(0,V,ev,'max'))
                ratio = helper(tt_delta);
                temp = [temp,ratio];
                
                temp_data = filttheta;
                [U,S,V,r,tspan,x,y,t] = hankle_analysis_v2(temp_data',dt, stackmax, rmax);
                V = V(:,1:r);
                ev = diag(S);
                ev = ev(1:r,1);
                ev = ev./sum(ev);
                t_delta = system_recon(sig_mean,V,ev,method);
                ratio = helper(t_delta);
                temp = [temp,ratio];
                
                temp_data = filtalpha;
                [U,S,V,r,tspan,x,y,t] = hankle_analysis_v2(temp_data',dt, stackmax, rmax);
                V = V(:,1:r);
                ev = diag(S);
                ev = ev(1:r,1);
                ev = ev./sum(ev);
                tt_alpha = system_recon(sig_mean,V,ev,method);
                ratio = helper(tt_alpha);
                temp = [temp,ratio];
                
                temp_data = filtbeta;
                [U,S,V,r,tspan,x,y,t] = hankle_analysis_v2(temp_data',dt, stackmax, rmax);
                V = V(:,1:r);
                ev = diag(S);
                ev = ev(1:r,1);
                ev = ev./sum(ev);
                tt_beta = system_recon(sig_mean,V,ev,method);
                ratio = helper(tt_beta);
                temp = [temp,ratio];
                
                temp_data = filtgamma;
                [U,S,V,r,tspan,x,y,t] = hankle_analysis_v2(temp_data',dt, stackmax, rmax);
                V = V(:,1:r);
                ev = diag(S);
                ev = ev(1:r,1);
                ev = ev./sum(ev);
                tt_gamma = system_recon(sig_mean,V,ev,method);
                ratio = helper(tt_gamma);
                temp = [temp,ratio];
                
                
            end
            res(window,:) = temp;
        end
        
        %         figure
        %         plot(mean(res(:,:,1),2))
        %         hold on
        %         plot(std(res(:,:,1),0,2))
        %         plot(mean(res(:,:,2),2))
        %         plot(std(res(:,:,2),0,2))
        %         legend('mean_mean','mean_std','std_mean','std_std','Location','northeastoutside')
        %         ylim=get(gca,'Ylim');
        %         xlabel('Time window')
        
        %         for s = 1:length(index_seizure)
        %             plot([index_seizure(s),index_seizure(s)],ylim,'m--')
        %         end
        %         filename = strrep(filename,'.mat','_attractor_');
        %         title([filename,'mean'])
        %
        %
        %         saveas(gcf,['./figures/',filename,num2str(stackmax),'_',num2str(rmax),'_',num2str(overlap*100),'.png'])
        %         close(gcf)
        
        filename = strrep(filename,'.mat','');
        mainname = ['./result_seperate_allbands_max_v4_',num2str(windows_second),'s_',num2str(stackmax),'_',num2str(rmax),'_',num2str(overlap*100),'_',num2str(seizure_percent*100),'/'];
        mkdir([mainname ,foldername])
        save([mainname,foldername,'/',filename,'.mat'],'res','label','windows_second','stackmax','overlap','seizure_percent','sessions')
        tline = fgetl(fid);
        
    end
end


function ratio = helper(tt)
%GETRECONFEA Summary of this function goes here
%   Detailed explanation goes here
ratio = [max(tt),min(tt),mean(tt),std(tt),mean(tt.^2),entropy(tt)];
% ratio = [ratio,modwtvar(modwt(tt,'db1',6),'db1')'];
end


