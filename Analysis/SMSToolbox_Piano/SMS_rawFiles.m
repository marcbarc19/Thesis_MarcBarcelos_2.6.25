function OUTfilenames = SMS_rawFiles(INfilenames) 



rawfilenum=1;
for curFile = 1:length(INfilenames)
    
    nameParts = strsplit(INfilenames(curFile).name,'_');
    
    if length(nameParts)==4
        
        OUTfilenames(rawfilenum) = INfilenames(curFile);
        rawfilenum = rawfilenum+1;
        
        
    end
end

