# Using Mendelian randomization to assess long-term pleiotropic effects of potential novel triglyceride-lowering medications

In the current study, by leveraging two large EHR-based biobanks, we examined genetic variants reported previously in MR studies from 5 TG-lowering target genes and tested their associations with clinical phenotypes.

**Background:** Drugs targeting triglyceride (TG)-associated genes have the potential to improve cardiovascular outcomes for patients with elevated TG levels. However, we know little regarding the potential additional benefits or deleterious effects of such targeting, particularly among individuals of African ancestry (AA). Mendelian randomization and PheWAS approaches offer the opportunity to examine such primary and secondary effects.

**Methods:** We examined 12 variants reported previously in Mendelian randomization studies from 5 genes that have been identified as TG-lowering targets (APOA5, LPL, APOC3, ANGPTL3, and ANGPTL4); for those variants associated with measured TG levels, we tested selected phenotypes, including lipid, cardiovascular, and other potential effects reported in previous studies, using PheWAS in separate cohorts of European ancestry (EA) patients and AA patients in BioVU. We also tested unspecified other phenotypes (i.e., without previously reported associations with TGs) for additional effects. We then replicated results in All of Us (AoU). As a secondary analysis, we tested the genetically predicted expression of these TG-lowering target genes for their association with the selected phenotypes in EA BioVU patients.

**Results:** Among BioVU EA patients (n=63,094), 11 previously reported SNPs were associated with measured TGs; of these, 9 SNPs were associated with lipid and cardiovascular phenotypes. Results were largely consistent in AoU EA participants (n=97,532). Among AA patients in BioVU (n=12,515) and AoU (n=31,710), results were more limited; only 6 of the 12 reported SNPs were associated with measured TGs in BioVU AA patients. While 4 of these validated 6 SNPs were associated with a lipid or cardiovascular phenotype in either BioVU or AoU, none were consistent across both cohorts. Additionally, we detected few secondary effects in either EA or AA BioVU patients, and none were replicated. In the secondary analysis assessing predicted gene expression, results were largely consistent with the primary analysis for EA BioVU patients.

**Conclusions:** These results suggest that beyond cardiovascular benefits there may be limited additional benefits, but few deleterious effects, from targeting known TG-associated genes for individuals of EA. However, we found limited information supporting the efficacy or safety of these targets for mitigating cardiovascular risk among AA individuals.

The included files execute the following within BioVU:
1) Run_Phewas_EA_cohort: conducts the PheWAS for validated SNPs in the EA BioVU cohort
2) Run_Phewas_AA_cohort: conducts the PheWAS for validated SNPs in the AA BioVU cohort
3) Run_Phewas_tissue: conducts the PheWAS for genetically predicted expression in selected tissues for 5 TG-associated genes
4) GBJ_analysis: conducts the Generalized Berk-Jones analysis for the TG-associated genes
