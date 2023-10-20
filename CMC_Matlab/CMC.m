
function rwu = CMC(tensor)

tensor = double(tensor);
nM = numel(size(tensor));
Margin = cell([nM 1]);
for iM = 1:nM
    iM_ = setdiff([nM:-1:1],iM);
    tmp = sum(tensor,iM_);
    Margin{iM} = tmp(:);
end

rwu = CMC_V0(Margin);



