prec_days.RData') # precs=[locs,98yrs,30d]
locs=cbind(lon,lat)

load(file='data/prec_all.RData') # lon,lat,prec

We analyzed two datasets in our article. 
Both are subsets of total precipitation rate (in m/s) files, which are publicly available and produced by the Community Earth System Model (CESM) Large Ensemble Project. 
Specifically, we considered log-transformed total precipitation rate in July 1 in 98 consecutive years starting in the year 402. 
For ease of working with the data and reproducing our results, we provide the two datasets in a GitHub repository accompanying our article.
