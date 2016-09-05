function reMakeMatfile(filePath)
% In some case, matfile() fails to access the saved variables, but load()
% can access them. This function make a new matfile with the same fields.

if strcmp(filePath(end-3:end),'.mat'),
    tmpFilePath = [filePath(1:end-3) '_tempSaveReMakeMatfile.mat'];
else
    tmpFilePath = [filePath '_tempSaveReMakeMatfile.mat'];
end
dataPtNew = matfile(tmpFilePath);
vars = load(filePath,'-mat');
names = fieldnames(vars);
for c = 1:length(names),
    dataPtNew.(names{c}) = vars.(names{c});
end
delete(filePath);
if strcmp(filePath(end-3:end),'.mat'),
    movefile(tmpFilePath,filePath);
else
    movefile(tmpFilePath,[filePath '.mat']);
end