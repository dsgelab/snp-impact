#!/usr/bin/env python
# coding: utf-8

import hail as hl
hl.init()

vcf_list = ['gs://finngen-production-library-red/finngen_R6/genotype_1.0/data/finngen_R6_chr%s.vcf.gz' % chrom for chrom in range(1, 24)]
mt = hl.import_vcf(vcf_list, force_bgz=True, reference_genome='GRCh38')

t = hl.import_table('gs://dsgelab-sakari/finemap_rsids_and_variant_ids_fixed_29072020.csv',
                    impute=True, delimiter=",", quote='"')
loci_to_extract = [hl.parse_locus('chr'+s, reference_genome='GRCh38') for s in t.locus_hg38.collect()]

mt_f = mt.filter_rows(hl.literal(loci_to_extract).contains(mt.locus))
mt_f = mt_f.key_rows_by('locus', 'alleles', 'rsid')
mt_f.GT.n_alt_alleles().export('gs://dsgelab-sakari/doc_jukarainen_SNPs_R6_25082020.tsv')
