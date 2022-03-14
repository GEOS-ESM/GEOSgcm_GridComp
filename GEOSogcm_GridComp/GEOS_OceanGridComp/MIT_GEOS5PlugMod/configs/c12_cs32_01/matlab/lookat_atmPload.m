% Code crashes because of GEOS/MITgcm land mask inconsistency when
% "pLoad(i,j,bi,bj)=p%PS(iSLo+i,jSLo+j)" in import_state_fill_mod.FOR
% pLoad is set to zero for time being.

% To test pLoad results, uncomment "pLoad(i,j,bi,bj)=p%PS(iSLo+i,jSLo+j)"
% in mitgcm_setup/code_split_driver/state/import/import_state_fill_mod.FOR
% Compile and run c12_cs32_01, then look at atmPload output files.

% Script was tested for compatibility with Octave 4.2.2.

% Make sure WorkingDir, MITgcm, and TEST directories are set correctly
WorkingDir='~/geos5/';
eval(['addpath ' WorkingDir 'MITgcm/utils/matlab/cs_grid/read_cs'])
eval(['cd ' WorkingDir 'TEST/scratch/mitocean_run'])

% Plot MITgcm land mask
hFacC=read_cs_bin('hFacC.data',1,'real*4',32);
figure(1), mycrossmap(hFacC)

% Plot raw atmospheric pressure received from GEOS
P=read_cs_bin('atmPload.0000072006.data',1,'real*4',32);
figure(2), mycrossmap(P)

% Atmospheric pressure with pressure over GEOS land mask
% set to mean atmospheric pressure over GEOS wet cells.
P2=P;
P2(find(P>1e14))=mean(P(find(P<1e14)));
figure(3), mycrossmap(P2)

% Atmospheric pressure with MITgcm land mask applied.
P3=P;
P3(find(hFacC==0))=mean(P(find(P<1e14)));
figure(4), mycrossmap(P3)
