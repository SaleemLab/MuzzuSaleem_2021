function waveforms = load_maxwaveforms(rfs, paradigm)

server = 'R:\viscog\';


for i=1:length(rfs)
    
    clusterID(i) = rfs{i}.KVclusterID;
    animal(i,:) = rfs{i}.animal;
    
end

[uanimals, ua, ub] = unique(animal,'rows');


for j=1:length(uanimals)
    
    
    %Look at xml files to decide which waveform file to use
    allXMLFiles = dir(fullfile(server,'EXPERIMENTS\UMA_DATA\',upper(uanimals(j,1:6)),uanimals(j,8:end),'*.xml'));
    allXMLFilesNames = {};
    [allXMLFilesNames{1:length(allXMLFiles)}] = deal(allXMLFiles.name);
    t = zeros(1, length(allXMLFilesNames));
    for tr = 1:length(allXMLFilesNames)
        t(tr) = ~isempty(strfind(lower(allXMLFilesNames{tr}),lower(paradigm)));
    end
    target = find(t==1);
    if ~isempty(target)
         targetFileName = regexprep(allXMLFilesNames{target},'[(.+)','');
         waveformfile = [targetFileName, '_Waveforms.mat'];   
    end
    
    Dat = load(fullfile(server,'EXPERIMENTS\UMA_DATA\',upper(uanimals(j,1:6)),uanimals(j,8:end),waveformfile));
    waveformDat{j} = Dat.waveformDat;
end    
    
for i=1:length(rfs)
    
    thiswaveformDat = waveformDat{ub(i)};
    icluster = find(cell2mat(thiswaveformDat.clusterIDs)==clusterID(i));
    
    thesewaveforms = thiswaveformDat.meanWaveform{icluster};
    
    %Get waveform from electrode site with largest amplitue
    [m im]= max(max(thesewaveforms'));
    waveforms(i,:) = thesewaveforms(im,:);
    
end

waveforms = -waveforms';