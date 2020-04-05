#./add_stat_shapes.py --filter diboson --prefix diboson_bbb WVchannel_datacard_BBB2.root WVchannel_datacard_BBB3.root"

./add_stat_shapes.py --filter DYJets_mWjj --prefix DYJets_bbb comb_XX.root comb_XX_tmp1.root
./add_stat_shapes.py --filter DYJets_mWV  --prefix DYJets_bbb comb_XX_tmp1.root comb_XX_tmp2.root

./add_stat_shapes.py --filter Top_mWjj --prefix Top_bbb comb_XX_tmp2.root comb_XX_tmp3.root
./add_stat_shapes.py --filter Top_mWV  --prefix Top_bbb comb_XX_tmp3.root comb_XX_tmp4.root

./add_stat_shapes.py --filter VBF_EWK_mWjj --prefix VBF_EWK_bbb comb_XX_tmp4.root comb_XX_tmp5.root
./add_stat_shapes.py --filter VBF_EWK_mWV  --prefix VBF_EWK_bbb comb_XX_tmp5.root comb_XX_tmp6.root

./add_stat_shapes.py --filter VBF_QCD_mWjj --prefix VBF_QCD_bbb comb_XX_tmp6.root comb_XX_tmp7.root
./add_stat_shapes.py --filter VBF_QCD_mWV  --prefix VBF_QCD_bbb comb_XX_tmp7.root comb_XX_tmp8.root

./add_stat_shapes.py --filter WJets_mWjj --prefix WJets_bbb comb_XX_tmp8.root comb_XX_tmp9.root
./add_stat_shapes.py --filter WJets_mWV  --prefix WJets_bbb comb_XX_tmp9.root comb_XX_bbb.root
