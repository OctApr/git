function testPair
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'codes'),fullfile(root,'codes','my_function'));
C=configLoad; setupParser(C);
fa=fullfile(root,'data','csi (1)一人','rx_2_260904_171616.csi');
fb=fullfile(root,'data','csi(2)二人+无人','rx_2_260904_180803.csi');
% Official import uses interpolation, matching normal drag-and-drop imports.
a=read_rxs_log(fa,80); b=read_rxs_log(fb,80);
ma=parseRXSBundle(a); mb=parseRXSBundle(b);
out=fullfile(root,'tmp','pair_tests');
R1=runPair(a,b,'show_figures',false,'path_res',out);
R2=runPair(ma,mb,'show_figures',false,'path_res',out);
for j=1:2
    assert(isequaln(R1(j).windows,R2(j).windows));
    assert(R1(j).audit.from_workspace && R2(j).audit.from_workspace);
    assert(R2(j).audit.corrupt==0 && R2(j).audit.frames==80);
    assert(isequal(R2(j).metadata.timestamp_ns,R1(j).metadata.timestamp_ns));
end
R3=runPair(fa,fb,'show_figures',false);
old=load(fullfile(root,'results','20260908_113643_206','results.mat'),'Results');
for j=1:2
    name=R3(j).name;
    for k=1:numel(old.Results)
        [~,oldName]=fileparts(old.Results{k}.file);
        if string(oldName)==name
            assert(isequaln(R3(j).windows.score,old.Results{k}.groups{1}.windows.score));
        end
    end
    assert(~R3(j).audit.from_workspace && R3(j).audit.interpolated_frames==0);
end
fprintf('PASS: pair paths, raw cells and merged workspace bundles.\n');
end
