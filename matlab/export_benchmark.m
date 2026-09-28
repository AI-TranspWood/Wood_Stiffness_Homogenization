% Export MATLAB reference results for the Python port (test_vs_matlab.py).
% Runs hom_TranspWood_model on the benchmark cases, stores every RVE level
% plus timing, and adds raw Hill tensors of fun_P_ellipsoid_aniso.

here = fileparts(mfilename('fullpath'));
stiff = fullfile(here,'..','..','02_stiffness');
addpath(stiff);

cases = {
 'birch_s1',       struct('wood','Birch','density_kg_m3',646,'moisture_percent',10,'crystallinity',0.58,'MFA_deg',10,'material_state',1)
 'birch_s2_HEMA',  struct('wood','Birch','density_kg_m3',646,'moisture_percent',10,'crystallinity',0.58,'MFA_deg',10,'material_state',2,'polymer','HEMA')
 'birch_s3_HEMA',  struct('wood','Birch','density_kg_m3',646,'moisture_percent',10,'crystallinity',0.58,'MFA_deg',10,'material_state',3,'polymer','HEMA')
 'birch_sweep',    struct('wood','Birch','MFA_deg',[0 5 10 20 30],'material_state',1)
 'birch_IF',       struct('wood','Birch','MFA_deg',10,'material_state',1,'interface_compliance',[0.1,0.05])
 'birch_s2_PMMA',  struct('wood','Birch','MFA_deg',10,'material_state',2,'polymer','PMMA')
 'spruce_s1',      struct('wood','Spruce','material_state',1)
 'spruce_s2_HEMA', struct('wood','Spruce','material_state',2,'polymer','HEMA')
 'spruce_sweep',   struct('wood','Spruce','MFA_deg',[0 10 20 30],'material_state',1)
 'pine_s1',        struct('wood','Pine','material_state',1)
 };

levels = {'pn','cel','cw','EWuc','LWuc','vesselwood','ringwood','ray','clearwood'};
bench = struct();
for ic = 1:size(cases,1)
    name = cases{ic,1};
    t = tic;
    r = hom_TranspWood_model(cases{ic,2});
    bench.(name).time = toc(t);
    bench.(name).MFA = r.MFA_deg;
    bench.(name).C_clearwood = r.C.clearwood;
    bench.(name).C_cellwall = r.C.cellwall;
    bench.(name).E = [r.E.R; r.E.T; r.E.L];
    bench.(name).G = [r.G.TL; r.G.RL; r.G.RT];
    bench.(name).nu = [r.nu.TR; r.nu.LR; r.nu.TL; r.nu.LT];
    bench.(name).M_cw = r.indentation.cellwall;
    for im = 1:numel(r.MFA_deg)
        for il = 1:numel(levels)
            if isfield(r.RVE{im},levels{il})
                bench.(name).(['L_',levels{il}])(:,:,im) = r.RVE{im}.(levels{il}).Chom;
            end
        end
    end
    fprintf('%-16s %8.2f s\n',name,bench.(name).time)
end

% Raw Hill tensors in real anisotropic matrices (EWuc and clearwood C0)
r = hom_TranspWood_model(cases{1,2});
hill.C0 = cat(3,r.RVE{1}.EWuc.C0,r.RVE{1}.clearwood.C0,r.RVE{1}.cw.Chom);
hill.shapes = [1 1; 1 0.01; 1 1e-3; 10 1e-20; 1 1e20; 3 0.2];
hill.P = zeros(6,6,size(hill.shapes,1),3);
t = tic;
for k = 1:3
    for s = 1:size(hill.shapes,1)
        hill.P(:,:,s,k) = fun_P_ellipsoid_aniso(hill.C0(:,:,k), ...
            hill.shapes(s,1),hill.shapes(s,2),8);
    end
end
hill.time_per_call = toc(t)/(3*size(hill.shapes,1));
fprintf('fun_P_ellipsoid_aniso: %.3f s per call\n',hill.time_per_call)

% Cross-check: hom_TranspWood_final (state 1) on a temporary copy that
% writes its result into tempdir instead of 02_stiffness (eval inside
% 02_stiffness so its relative addpath works).
src = fileread(fullfile(stiff,'hom_TranspWood_final.m'));
out = fullfile(tempdir,'res_final_bench.mat');
src = strrep(src,'save(output_file,''results'');', ...
    ['save(''',out,''',''results'');']);
src = strrep(src,'clear;','');
old = cd(stiff);
t = tic; eval(src); final_time = toc(t);
cd(old);
S = load(out);
final.C_clearwood = S.results(1).RVE.clearwood.Chom;
final.time = final_time;

save(fullfile(here,'..','bench_matlab.mat'),'bench','hill','final','-v7');
fprintf('saved bench_matlab.mat\n')
