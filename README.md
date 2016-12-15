# heavyhitters

Scripts to download and analyze heavy hitter concentration/stability statistics on CAIDA passive-2016 anonymized Internet trace data

Plots:
core_correlation_(minute/hour).png: density of intersection of core list in time window t with core list in time window 0
cover_corerlation_(minute/hour).png: percentage of sum of all edge weights in graph for window t that are covered by core list in time window 0
cover_densities_(minute/hour).png: (95% cover size)/(total # distinct IPs) in one time window