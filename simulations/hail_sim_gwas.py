#!/usr/bin/env python
# coding: utf-8
​
​
​
import hail as hl
import random
hl.init()
​
​
# t = hl.read_matrix_table('gs://mattia/mattia-simulations/simEUR350_correlated.mt')
# t = t.drop(t.beta,t.y,t.sex,t.y_no_noise,t.y,t.ldscsim)
# pruned = hl.ld_prune(t.GT, r2=0.1, bp_window_size=1000000)
# t_pruned = t.filter_rows(hl.is_defined(pruned[t.row_key]))
# t_pruned.write('gs://mattia/output_sakari_sim/simEUR350_correlated_pruned.mt')
​
​
def gwas(y, x, cov):
    g = hl.linear_regression_rows(y=y,
                                  x=x,
                                  covariates=cov,
                                  pass_through=['rsid'])
    return g
​
def gwas_logistic(y, x, cov):
    g = hl.logistic_regression_rows(test='wald',
                                  y=y,
                                  x=x,
                                  covariates=cov,
                                  pass_through=['rsid'])
    return g
​
​
t = hl.read_matrix_table('gs://mattia/output_sakari_sim/simEUR350_correlated_pruned.mt')
​
for i in range(0,72):
    print(i)
    random.seed(123+i)
    h2=random.randrange(10,60,1)/100
    pi=random.randrange(1,20,1)/100
    k=random.randrange(1,20,1)/100
    random.seed(123+i)
    sim=hl.experimental.ldscsim.simulate_phenotypes(t,t.GT,h2=h2,pi=pi)
    sim=hl.experimental.ldscsim.binarize(sim,y=sim.y,K=k)
​
    to_exp = sim.annotate_rows(AF=sim.variant_qc.AF[0]).rows().key_by().select('rsid','AF','beta')
    to_exp.export('gs://mattia/output_sakari_sim/rsid_AF_beta_gold_pheno' + str(i) + '_h2' +str(h2) + '_pi' + str(pi) + '_k' + str(k) + '.tsv.bgz')
    
    gw_results = gwas(sim.y, sim.GT.n_alt_alleles(), [1.0,sim.sample_info.PC1,sim.sample_info.PC2,sim.sample_info.PC3,sim.sample_info.PC4,sim.sample_info.PC5])
    to_exp = gw_results.key_by()
    to_exp.select('rsid','beta','standard_error','p_value').export('gs://mattia/output_sakari_sim/rsid_beta_se_gwas_pheno' + str(i) + '.tsv.bgz')
​
    gw_results_logi = gwas_logistic(sim.y_binarized, sim.GT.n_alt_alleles(), [1.0,sim.sample_info.PC1,sim.sample_info.PC2,sim.sample_info.PC3,sim.sample_info.PC4,sim.sample_info.PC5])
    to_exp_logi = gw_results_logi.key_by()
    to_exp_logi.select('rsid','beta','standard_error','p_value').export('gs://mattia/output_sakari_sim/rsid_beta_se_gwas_pheno' + str(i) + '_logistic.tsv.bgz')