for i = 1:numfiles
    filename = ['relativeZ_' num2str(i)];
    info = geotiffinfo(['container_' int2str(l) '_mola128.tif']);
    geotiffwrite(filename,diffM(i).d,r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
end
