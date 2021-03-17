#!/usr/bin/env python
# coding: utf-8

import hail as hl
hl.init()

vcf_list = ['gs://finngen-production-library-red/finngen_R7/genotype_2.0/data/finngen_R7_chr%s.vcf.gz' % chrom for chrom in range(1, 24)]
mt = hl.import_vcf(vcf_list, force_bgz=True, reference_genome='GRCh38')

mt.describe()

t = hl.import_table('gs://dsgelab-sakari/finemap_variants_18102020.tsv',
                    impute=True, delimiter="\t", quote='"')
loci_to_extract = [hl.parse_locus('chr'+s, reference_genome='GRCh38') for s in t.locus_hg38.collect()]

mt_f = mt.filter_rows(hl.literal(loci_to_extract).contains(mt.locus))
mt_f = mt_f.key_rows_by('locus', 'alleles', 'rsid')

mt_f.GT.n_alt_alleles().export('gs://dsgelab-sakari/snp_17032021/GT_17032021.tsv')
mt_f.DS.export('gs://dsgelab-sakari/snp_17032021/DS_17032021.tsv')
mt_f.rows().info.INFO.export('gs://dsgelab-sakari/snp_17032021/INFO_17032021.tsv')
