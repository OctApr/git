function [Groups, Audit] = fileLoad(filepath,Config)
%FILELOAD 按源 MAC 和物理配置分组；输出维度 [tx,rx,pkt,sc]。
% 逐帧边界检查；官方解码器关闭插值，避免人为引入子载波相关性。
if nargin<2, Config=configLoad; end
fromWorkspace=~(ischar(filepath) || (isstring(filepath)&&isscalar(filepath)));
if fromWorkspace
    frames=workspaceFrames(filepath); frameIndex=0;
else
    fid=fopen(filepath,'rb','ieee-le');
    assert(fid>=0,'RapidPD:FileOpen','无法打开 %s',filepath);
    cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
end
Audit=struct('frames',0,'decoded',0,'corrupt',0,'truncated',0, ...
    'missing_timestamp',0,'invalid_csi',0,'duplicate_time',0,'filtered',0, ...
    'from_workspace',fromWorkspace,'interpolated_frames',0, ...
    'messages',strings(0,1));
Groups=struct([]); keys=strings(0,1); firstTime=[];
while true
    if fromWorkspace
        frameIndex=frameIndex+1;
        if frameIndex>numel(frames), break; end
    else
    head=fread(fid,4,'*uint8');
    if isempty(head), break; end
    if numel(head)~=4, Audit.truncated=Audit.truncated+1; break; end
    n=double(typecast(head,'uint32'));
    if n<7 || n>64*1024*1024
        Audit.corrupt=Audit.corrupt+1;
        Audit.messages(end+1)="Invalid frame length; stopped at byte "+string(ftell(fid)-4);
        break;
    end
    body=fread(fid,n,'*uint8');
    if numel(body)~=n, Audit.truncated=Audit.truncated+1; break; end
    end
    Audit.frames=Audit.frames+1;
    try
        if fromWorkspace, r=frames{frameIndex}; else, r=RXSParser([head;body],false); end
        if isempty(r) || ~isscalar(r) || ~isfield(r,'CSI') || isempty(r.CSI)
            error('RapidPD:Frame','Missing/split CSI frame');
        end
        Audit.decoded=Audit.decoded+1;
        c=r.CSI;
        if isfield(c,'PhaseSlope') && ~isempty(c.PhaseSlope)
            Audit.interpolated_frames=Audit.interpolated_frames+1;
        end
        % MEX exposes counts as uint8/uint16: cast before products to avoid saturation.
        numericFields={'NumTx','NumRx','NumCSI','NumESS','PacketFormat','CBW', ...
            'CarrierFreq','FirmwareVersion','ANTSEL'};
        for nf=1:numel(numericFields), c.(numericFields{nf})=double(c.(numericFields{nf})); end
        mac=string(lower(join(compose('%02x',r.StandardHeader.Addr2),':')));
        if (~isempty(Config.packet_formats) && ~ismember(c.PacketFormat,Config.packet_formats)) || ...
                (strlength(string(Config.source_mac))>0 && mac~=lower(string(Config.source_mac)))
            Audit.filtered=Audit.filtered+1; continue;
        end
        if ~isfield(r,'RxSBasic') || ~isfield(r.RxSBasic,'SystemTime') || r.RxSBasic.SystemTime==0
            Audit.missing_timestamp=Audit.missing_timestamp+1; continue;
        end
        stamp=uint64(r.RxSBasic.SystemTime);
        if isempty(firstTime), firstTime=stamp; end
        x=c.CSI;
        sc=double(c.SubcarrierIndex(:));
        if c.NumCSI~=1 || c.NumESS~=0 || numel(x)~=numel(sc)*c.NumTx*c.NumRx || ...
                any(~isfinite(x(:))) || any(sum(abs(reshape(x,numel(sc),[])),1)==0)
            Audit.invalid_csi=Audit.invalid_csi+1; continue;
        end
        x=reshape(x,numel(sc),c.NumTx,c.NumRx);
        [sc,order]=sort(sc); x=x(order,:,:);
        % Include firmware, center, offset, antenna selection and exact tone set.
        key=sprintf('%s_f%d_bw%d_fc%.0f_tx%d_rx%d_fw%d_ant%d_sc%s',char(mac), ...
            c.PacketFormat,c.CBW,c.CarrierFreq,c.NumTx,c.NumRx,c.FirmwareVersion,c.ANTSEL, ...
            sprintf('%d,',sc));
        gi=find(keys==string(key),1);
        if isempty(gi)
            gi=numel(keys)+1;
            g=struct('source_mac',mac,'packet_format',c.PacketFormat,'bandwidth_mhz',c.CBW, ...
                'carrier_hz',c.CarrierFreq,'firmware',c.FirmwareVersion, ...
                'ntx',c.NumTx,'nrx',c.NumRx,'subcarrier_index',sc, ...
                'packet_csi',{{}},'timestamp_ns',uint64([]),'frame_number',[], ...
                'time_seconds',[],'csi_raw',[],'median_rate_hz',NaN);
            if isempty(Groups), Groups=g; else, Groups(gi)=g; end %#ok<AGROW>
            keys(gi)=string(key);
        end
        Groups(gi).packet_csi{end+1}=single(reshape(x,numel(sc),c.NumTx,c.NumRx));
        Groups(gi).timestamp_ns(end+1)=stamp;
        Groups(gi).frame_number(end+1)=Audit.frames;
    catch ME
        Audit.corrupt=Audit.corrupt+1;
        if numel(Audit.messages)<10, Audit.messages(end+1)=string(ME.message); end
    end
end
for gi=1:numel(Groups)
    g=Groups(gi);
    [st,ix]=sort(g.timestamp_ns);
    [st,uniqueIndex]=unique(st,'stable'); ix=ix(uniqueIndex);
    Audit.duplicate_time=Audit.duplicate_time+numel(g.timestamp_ns)-numel(st);
    % uint64 subtraction BEFORE double avoids loss of nanosecond precision.
    if any(st<firstTime)
        t=zeros(size(st));
        hi=st>=firstTime;
        t(hi)=double(st(hi)-firstTime)/1e9;
        t(~hi)=-double(firstTime-st(~hi))/1e9;
    else
        t=double(st-firstTime)/1e9;
    end
    a=cat(4,g.packet_csi{ix});
    Groups(gi).csi_raw=permute(a,[2 3 4 1]);
    Groups(gi).timestamp_ns=st(:);
    Groups(gi).frame_number=g.frame_number(ix).';
    Groups(gi).time_seconds=t(:);
    if numel(t)>1, Groups(gi).median_rate_hz=1/median(diff(t)); end
    Groups(gi).packet_csi={};
end
end
