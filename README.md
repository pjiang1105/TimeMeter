# TimeMeter
TimeMeter is a statistical tool and R package to assess temporal gene expression pattern similarity, and identify differential progressing genes. TimeMeter uses the dynamic time warping (DTW) algorithm to align two temporal sequences. TimeMeter first post-processes DTW alignment by truncating certain start or end points based on alignment patterns. This will result in a truncated alignment which represents the time frames from each time series that are comparable. Then it calculates four measurements that jointly assess gene pair temporal similarity: percentage of alignment for (1) query and for (2) reference, respectively; (3) aligned gene expression correlation; (4) likelihood of alignment arising by chance. These four novel metrics in TimeMeter give temporal similarity assessments on different aspects, and the joint requirement of these four metrics will give a high confident assessment for the temporal pattern similarity. For gene pairs with similar temporal patterns, TimeMeter partitions the temporal associations into separate segments via piecewise regression. The differential progression between gene pairs is calculated by aggregation of progression difference in each segment.

# Reference
Jiang, Peng, et al. TimeMeter assesses temporal gene expression similarity and identifies differentially progressing genes. Nucleic Acids Research (2020).

# Contact
Peng Jiang, Ph.D </br>
Assistant Professor </br>
Center for Gene Regulation in Health and Disease (GRHD) </br>
Center for Applied Data Analysis and Modeling (ADAM) </br>
Department of Biological, Geological and Environmental Sciences (BGES) </br>
Cleveland State University, 2121 Euclid Ave, Cleveland, OH 44115 </br>

Lab Website: https://sites.google.com/view/jiang-lab/
