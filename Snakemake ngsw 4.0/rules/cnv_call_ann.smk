#there is the cohort mode and the case mode

#The cohort mode simultaneously generates a cohort model and calls CNVs for the cohort samples.

#The case mode analyzes a single sample against an already constructed cohort model.
#The same workflow steps apply to both targeted exome and whole genome sequencing (WGS) data. 

#The workflow is able to call both rare and common events and intelligently handles allosomal ploidies, i.e. cohorts of mixed male and female samples.

# -ProcessInterval
# -AnnotateIntervals + CollectReadCounts into FilterIntervals
# -DetermineGermlineContigPloidy
# -GermlineCNVCaller
# -PostProcessGermlineCNVCalls