function setupParser(Config)
% 使用官方 MEX；不修改全局 ParserPreference，不自动编译或下载。
if ~isempty(Config.parser_path)
    addpath(Config.parser_path,fullfile(Config.parser_path,'parser_basic'));
end
if exist('RXSParser','file')~=3
    local=fullfile(getenv('USERPROFILE'),'PicoScenes-MATLAB-Toolbox-Core');
    if isfolder(local), addpath(local,fullfile(local,'parser_basic')); end
end
assert(exist('RXSParser','file')==3,'RapidPD:ParserMissing', ...
    ['找不到官方 RXSParser MEX。请安装 PicoScenes-MATLAB-Toolbox-Core，' ...
     '然后使用 main(''parser_path'',''工具箱目录'')。']);
end
