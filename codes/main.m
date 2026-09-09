function [Summary, Results] = main(varargin)
%MAIN 从任意目录调用：addpath('codes'); [summary,results]=main;
% 默认只输出分数；main('threshold',0.43) 为未校准的演示判决。
codeDir=fileparts(mfilename('fullpath')); addpath(fullfile(codeDir,'my_function'));
Config=configLoad(varargin{:}); setupParser(Config);
files=dir(fullfile(Config.path_data,'**',Config.file_pattern));
assert(~isempty(files),'RapidPD:NoFiles','%s 内未找到 .csi 文件。',Config.path_data);
runDir=fullfile(Config.path_res,char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
mkdir(runDir); Results=cell(numel(files),1); summaryRows=cell(0,13);
for fi=1:numel(files)
    filepath=fullfile(files(fi).folder,files(fi).name);
    relative=erase(string(filepath),string(Config.path_data)+filesep);
    fprintf('[%d/%d] %s\n',fi,numel(files),relative);
    [groups,audit]=fileLoad(filepath,Config);
    fileResult=struct('file',relative,'audit',audit,'groups',{{}});
    if isempty(groups)
        summaryRows(end+1,:)={relative,0,"",NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,"no_valid_frames"}; %#ok<AGROW>
    end
    for gi=1:numel(groups)
        g=groups(gi); wins=winSplit(g.time_seconds,Config); n=numel(wins);
        profile=string(sprintf('%s_fmt%d_bw%d_fc%.0f_tx%d_rx%d_sc%d', ...
            g.source_mac,g.packet_format,g.bandwidth_mhz,g.carrier_hz,g.ntx,g.nrx, ...
            numel(g.subcarrier_index)));
        groupConfig=Config;
        if strlength(string(Config.threshold_profile))>0 && profile~=string(Config.threshold_profile)
            groupConfig.threshold=NaN;
        end
        start=nan(n,1); stop=start; packets=zeros(n,1); maxGap=start;
        score=start; scoreMean=start; scoreSum=start; streamSum=start; raw=start;
        streams=nan(n,g.ntx*g.nrx); valid=false(n,1);
        for wi=1:n
            w=wins(wi); start(wi)=w.start; stop(wi)=w.stop;
            packets(wi)=numel(w.indices); maxGap(wi)=w.max_gap; valid(wi)=w.valid;
            if ~w.valid, continue; end
            d=winProcess(g.csi_raw(:,:,w.indices,:),groupConfig);
            score(wi)=d.score; scoreMean(wi)=d.score_mean; scoreSum(wi)=d.score_sum;
            streamSum(wi)=d.score_stream_sum; raw(wi)=d.present;
            streams(wi,:)=d.stream_statistics;
        end
        smooth=labelFilt(raw,Config.smooth_windows);
        windows=table(start,stop,packets,maxGap,valid,score,scoreMean,scoreSum,streamSum,raw,smooth);
        windows.profile=repmat(profile,n,1);
        for si=1:size(streams,2), windows.(sprintf('stream_%d',si))=streams(:,si); end
        metadata=rmfield(g,{'csi_raw','packet_csi'});
        fileResult.groups{gi}=struct('metadata',metadata,'windows',windows);
        name=sprintf('%02d_%s_group%02d',fi,erase(files(fi).name,'.csi'),gi);
        writetable(windows,fullfile(runDir,[name '.csv']));
        if Config.make_plots && n>0
            resDisp(windows,g,groupConfig,fullfile(runDir,[name '.png']));
        end
        status="ok";
        if ~any(valid), status="insufficient_data";
        elseif ~isfinite(groupConfig.threshold), status="scores_only_uncalibrated";
        else, status="threshold_supplied"; end
        meanScore=mean(score,'omitnan');
        presentFraction=mean(smooth,'omitnan');
        summaryRows(end+1,:)={relative,gi,g.source_mac,g.packet_format,g.bandwidth_mhz, ...
            size(g.csi_raw,3),numel(g.subcarrier_index),sum(valid),g.median_rate_hz, ...
            meanScore,presentFraction,audit.corrupt+audit.truncated,status}; %#ok<AGROW>
    end
    Results{fi}=fileResult;
end
Summary=cell2table(summaryRows,'VariableNames',{'file','group','source_mac','packet_format', ...
    'bandwidth_mhz','packets','subcarriers','valid_windows','median_rate_hz', ...
    'mean_score','presence_fraction','file_corrupt_frames','status'});
writetable(Summary,fullfile(runDir,'summary.csv'));
save(fullfile(runDir,'results.mat'),'Config','Summary','Results','-v7.3');
fid=fopen(fullfile(runDir,'config.json'),'w');
if fid>=0, fprintf(fid,'%s',jsonencode(Config,PrettyPrint=true)); fclose(fid); end
fprintf('完成：%s\n',runDir);
if isnan(Config.threshold)
    fprintf('尚未设置校准阈值：raw/smooth 为 NaN，不代表无人。\n');
end
fprintf('文件 %d，汇总行 %d，有效窗口 %d。\n',numel(files),height(Summary),sum(Summary.valid_windows));
disp(Summary(:,{'file','valid_windows','mean_score','status'}));
if Config.make_plots
    plotComparison(fullfile(runDir,'results.mat'),'','',Config.show_figures);
end
end
