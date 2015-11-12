for i = 1:numfiles
    info = geotiffinfo(['container_' int2str(i) '_mola128.tif']);
    
    geotiffwrite(['image_output/RadMC_C' num2str(i) '_M1.tif'],double(mat(i).mpolarr(1).mp(:,:,1)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(['image_output/AngMC_C' num2str(i) '_M1.tif'],double(mat(i).mpolarr(1).mp(:,:,2)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(['image_output/RelZ_C' num2str(i) '_M1.tif'],double(mat(i).mpolarr(1).mp(:,:,3)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(['image_output/FixZ_C' num2str(i) '_M1.tif'],double(mat(i).mpolara(1).mp(:,:,3)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    
    geotiffwrite(['image_output/AngCW_C' num2str(i) '_M1.tif'],double(mat(i).mpolarr(1).mp(:,:,1)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(['image_output/RadCW_C' num2str(i) '_M1.tif'],double(mat(i).mpolarr(1).mp(:,:,1)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    
    for nmic = 2:10
        if length(m(i).m)>=nmic
            geotiffwrite(['image_output/RadMC_C' num2str(i) '_M' num2str(nmic) '.tif'],double(mat(i).mpolarr(nmic).mp(:,:,1)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
            geotiffwrite(['image_output/AngMC_C' num2str(i) '_M' num2str(nmic) '.tif'],double(mat(i).mpolarr(nmic).mp(:,:,2)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
            geotiffwrite(['image_output/RelZ_C' num2str(i) '_M' num2str(nmic) '.tif'],double(mat(i).mpolarr(nmic).mp(:,:,3)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
            geotiffwrite(['image_output/FixZ_C' num2str(i) '_M' num2str(nmic) '.tif'],double(mat(i).mpolara(nmic).mp(:,:,3)),r(i).r,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
        end
    end
end