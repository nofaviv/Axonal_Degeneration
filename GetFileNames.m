function FoundFileNames=GetFileNames(SearchDir,FileType)

FoundFileNames=char([]);
if(isempty(SearchDir)||isempty(FileType))
    return;
end;

FileNames=dir(SearchDir);
for n=1:size(FileNames,1)    
    FullFileName=[SearchDir '\' FileNames(n).name];    
    if(FileNames(n).isdir)
        if(~strcmp(FileNames(n).name,'.') && ~strcmp(FileNames(n).name,'..'))
            NextLevelFileNames=GetFileNames(FullFileName,FileType);
            if(~isempty(NextLevelFileNames))
                FoundFileNames((size(FoundFileNames,1)+1):(size(FoundFileNames,1)+size(NextLevelFileNames,1)),1:size(NextLevelFileNames,2))=NextLevelFileNames;            
            end;
        end;
    elseif(~isempty(strfind(FileNames(n).name,['.' FileType])))        
        FoundFileNames(size(FoundFileNames,1)+1,1:size(FullFileName,2))=FullFileName;        
    end;
end;
