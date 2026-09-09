function Summary = auditRapidPD
%AUDITRAPIDPD 对指定的一人运动与无人数据复核第二轮 FACF 分数。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'codes','my_function'));
C=configLoad; setupParser(C);
files=["csi (1)一人/rx_2_260904_171616.csi", ...
       "csi(2)二人+无人/rx_2_260904_180803.csi"];
labels=["one_person_motion";"empty"];
rows=cell(numel(files),6); windowTables=cell(numel(files),1);
for fi=1:numel(files)
    path=fullfile(root,'data',files(fi)); [groups,audit]=fileLoad(path,C);
    assert(numel(groups)==1 && audit.corrupt==0);
    g=groups(1); wins=winSplit(g.time_seconds,C); score=nan(numel(wins),1);
    for wi=1:numel(wins)
        score(wi)=winProcess(g.csi_raw(:,:,wins(wi).indices,:),C).score;
    end
    W=table([wins.start].',[wins.stop].',repmat(20,numel(wins),1),score, ...
        'VariableNames',{'start_s','stop_s','packets','score'});
    W.file=repmat(files(fi),height(W),1); W.label=repmat(labels(fi),height(W),1);
    windowTables{fi}=W;
    rows(fi,:)={files(fi),labels(fi),size(g.csi_raw,3),numel(wins),mean(score),std(score)};
end
Summary=cell2table(rows,'VariableNames', ...
    {'file','label','packets','windows','mean_score','std_score'});
out=fullfile(root,'results','audit_round2'); if ~isfolder(out), mkdir(out); end
writetable(Summary,fullfile(out,'summary.csv'));
writetable(vertcat(windowTables{:}),fullfile(out,'windows.csv'));
disp(Summary);
end
