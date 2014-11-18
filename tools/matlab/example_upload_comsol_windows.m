mphstart();
m=mphload_enhanced([pwd,'\graphene_eigen_3d.mph']);
m.sol('sol1').feature('e1').active(false);
m.sol('sol1').feature('v1').set('scalemethod','none');

id_extract=tic();
mtrx=mphmatrix(m,'sol1','Out',{'Kc','Ec','Dc'});
fprintf('Extracting took %f seconds.\n',toc(id_extract));

id_write=tic();
PetscBinaryWrite('./comsol_gr3d/K.dat',sparse(mtrx.Kc),'ispetsccomplex',true);
PetscBinaryWrite('./comsol_gr3d/E.dat',sparse(mtrx.Ec),'ispetsccomplex',true);
PetscBinaryWrite('./comsol_gr3d/D.dat',sparse(mtrx.Dc),'ispetsccomplex',true);
fprintf('Writing data took %f seconds.\n',toc(id_write));

id_upload=tic();
PSCP_CMD_BASE='pscp -batch -q -r';
PSCP_CMD_KEY='-i C:\Users\Ian_2\Dropbox\secure\ian-tacc.ppk';
PSCP_CMD_SRC='.\comsol_gr3d';
PSCP_CMD_DEST='iwill@stampede.tacc.utexas.edu:/work/03165/iwill/';
status=system([PSCP_CMD_BASE,' ',PSCP_CMD_KEY,' ',PSCP_CMD_SRC,' ',PSCP_CMD_DEST]);
if status==0
    fprintf('Successfully uploaded to TACC in %f seconds.\n',toc(id_upload));
else
    fprintf('Failure in uploading to TACC.\n');
end

fprintf('Total run time was %f seconds.\n',toc(id_write));

