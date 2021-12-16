function generate_slurm_script()
str = fileread('X:\Lab\Labmembers\Urs Boehm\scripts\SLURM\PCA_ICA_array_template.slurm');
fid = fopen('PCA_ICA_array.slurm', 'w');
files = dir('*_rois.mat');
nfiles = length(files);
currPath = pwd;
currPath = strrep(currPath,'\','/');
currPath = strrep(currPath,'X:','/n/cohen_lab');
str = strrep(str,'datapath',currPath);
fprintf(fid,str,nfiles);
fclose(fid);
end