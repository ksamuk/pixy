from itertools import combinations
from collections import Counter
from scipy import special

#For the given region: return average pi, # of differences, # of comparisons, and # missing.
# this function loops over every site in a region passed to it
def tallyRegion(gt_region):
        total_diffs = 0
        total_comps = 0
        total_missing = 0
        for site in gt_region:
            vec = site.flatten()
            #now we have an individual site as a numpy.ndarray, pass it to the comparison function
            site_diffs, site_comps, missing = compareGTs(vec)
            total_diffs += site_diffs
            total_comps += site_comps
            total_missing += missing
        if total_comps > 0:
            avg_pi = total_diffs/total_comps
        else:
            avg_pi = 0
        return(avg_pi, total_diffs, total_comps, total_missing)

#For the given region: return average dxy, # of differences, # of comparisons, and # missing.
# this function loops over every site in a region passed to it
def dxyTallyRegion(pop1_gt_region, pop2_gt_region):
	total_diffs = 0
	total_comps = 0
	total_missing = 0
	for x in range(0,len(pop1_gt_region)):
		site1 = pop1_gt_region[x]
		site2 = pop2_gt_region[x]
		vec1 = site1.flatten()
		vec2 = site2.flatten()
		#now we have an individual site as 2 numpy.ndarrays, pass them to the comparison function
		site_diffs, site_comps, missing = dxyCompareGTs(vec1, vec2)
		total_diffs += site_diffs
		total_comps += site_comps
		total_missing += missing
	if total_comps > 0:
		avg_pi = total_diffs/total_comps
	else:
		avg_pi = 0
	return(avg_pi, total_diffs, total_comps, total_missing)

#Return the number of differences, the number of comparisons, and missing data count.
def compareGTs(vec): #for pi
	c = Counter(vec)
	diffs = c[1]*c[0]
	gts = c[1]+c[0]
	missing = (len(vec))-gts  #anything that's not 1 or 0 is ignored and counted as missing
	comps = int(special.comb(gts,2))
	return(diffs,comps,missing)

def dxyCompareGTs(vec1, vec2): #for dxy
	c1 = Counter(vec1)
	c2 = Counter(vec2)
	gt1zeros = c1[0]
	gt1ones = c1[1]
	gts1 = c1[1]+c1[0]
	gt2zeros = c2[0]
	gt2ones = c2[1]
	gts2 = c2[1]+c2[0]
	missing = (len(vec1)+len(vec2))-(gts1+gts2)  #anything that's not 1 or 0 is ignored and counted as missing  
	diffs = (gt1zeros*gt2ones)+(gt1ones*gt2zeros)
	comps = gts1*gts2
	return(diffs,comps,missing)