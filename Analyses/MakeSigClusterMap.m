function MakeSigClusterMap(cfg)
% function MakeSigClusterMap(cfg)
%
% cfg = [];
% cfg.root      = fullfile(root,'Results');
% cfg.outDir    = 'GroupResults\FirstLevel\MT_BehModRegressors\Viv&RJConj';
% cfg.prefix    = 'RJ_Viv';
% cfg.tmap      = 4;

outDir    = fullfile(cfg.root,cfg.outDir);
[V,tvals] = read_nii(fullfile(outDir,sprintf('spmT_000%d.nii',cfg.tmap)));

cluster_files = str2fullfile(outDir,[cfg.prefix, '*.nii']);
if ~iscell(cluster_files); tmp = cluster_files; clear cluster_files; cluster_files{1} = tmp; end
nCs = length(cluster_files); clusters = zeros(V.dim);
for c = 1:nCs; [~,tmp] = read_nii(cluster_files{c}); 
    clusters = clusters + tmp; clear tmp; end

sigT = zeros(V.dim); sigT(clusters>0) = tvals(clusters>0);
write_nii(V,sigT,fullfile(outDir,sprintf('sigT_000%d.nii',cfg.tmap)));