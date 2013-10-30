mkdir lib/drupke
cp ../spectra/ebv_ccm.pro lib/drupke
cp ../spectra/extcurve_ccm.pro lib/drupke
cp ../io/readcol300.pro lib/drupke
cp ../io/readspec.pro lib/drupke
cp ../misc/cmset_op.pro lib/drupke
cp ../science/gaussarea.pro lib/drupke
cp ../science/gaussflux.pro lib/drupke
cp ../misc/jjadd_tag.pro lib/drupke
tar zcvf uhspecfit_lris.tgz README_lris example* common lib lris stellar_models/gonzalezdelgado/SSPGeneva_z020.sav stellar_models/gonzalezdelgado/gd_template.pro stellar_models/gonzalezdelgado/README
rm -rf lib/drupke
