# DifFracTion
## Background 

Hi-C have significantly enhanced our understanding of the 3D organization of the genome. By mapping pairwise interactions between distal genomic regions, Hi-C enables the discovery of hierarchical 3D organization of the genome.

However, Hi-C contains many DNA sequence and technology driven biases, that prevent effective comparison of chromatin interactions between different conditions. 

Technology and sequence biases, such as the cutting frequency of the chosen restriction enzyme, cross-linking efficiency, and GC content, can be alleviated within a single Hi-C data set using matrix balancing algorithms. The core idea of matrix balancing is straightforward: if most genomic loci have roughly equal visibility, then differences in row or column sums of the contact matrix are more likely to reflect technical artifacts than biology. Balancing rescales the matrix, so that each row and column contributes equally, improving within-dataset normalization.

However, while matrix balancing corrects biases inside one dataset, it does not remove differences that exists between datasets generated under different conditions. Factors such as sequencing depth, library complexity, or global chromatin reorganization can cause systematic shifts across experiments. Therefore, when comparing Hi-C maps across cell types or conditions, additional strategies are needed to harmonize datasets. 

HiC experiments have shown that the contact probability between two loci $\mathbf{P_c}$, separated  separated by a genomic distance $\mathbf{s}$, follows a power law 

$$\mathbf{P_c(s) \sim 1/s^{\alpha}}$$
where the exponent ${\alpha}$ is a **scaling exponent** that describes the rate of distance decay.

Aiden *et al.* described this decay by averaging across all chromosomes in a human cell line over genomic distances between 0.5 Mb and 7 Mb, reporting an exponent of approximately $\alpha = 1.08$. Some chromosomes showed deviations from this average.

Deviations from the average scaling exponent may reflect genuine biological differences, such as altered chromatin organization, and should be preserved. However, technical or dataset specific biases can still distort the absolute counts that underlie the decay curve. 

To address this, we have developed DifFracTion, an algorithm that iteratively rescales Hi-C contact counts so that datasets with different technical biases can be placed on a comparable scale and compared. The central idea is to harmonize contact counts globally while leaving the scaling exponent $\alpha $ intact, ensuring that differences in $\alpha$ across conditions remain interpretable as biological variation rather than artifacts. 



## Iterative Rescaling of Distance Decay 
After Knight–Ruiz (KR) matrix balancing is applied to each Hi-C map, we implement a normalization step designed to harmonize contact counts across datasets while preserving their intrinsic distance-decay scaling. 
First, we subset genomic distances within the range known to capture the global decay trend (0.5-7Mb). For each dataset, we computed mean contact counts as a function of genomic distance. A 'reference/dataset A' is selected and fixed, while the remaining (dataset B) is iteratively rescaled to better match the reference.
At each iteration DifFracTion:
1) Computes distance specific ratios $\mathbf{r}$ : for every genomic distance $\mathbf{s}$, we calculate the ratio $\mathbf{r_s}$ between the dataset A and dataset B.
2) Weight by distance: To prioritize short-range interactions (which dominate Hi-C signal), ratios are down-weighted with a power-law kernel $$\mathbf{w(s) \propto s^{-1.08}}$$. We normalize weights to the range [0,1] to keep ratios intact while controlling their influence ($\mathbf{q = w_s * r_s}$).

3) Derive a correction factor: To avoid overly aggressive updates, we apply a smoothing parameter $\mathbf{\eta}$ , so that the effective weighted factor is:
 $\mathbf{f = 1 + \eta * (\text{q} - 1),}$
where $\mathbf{\eta}=0.15$ represents the fraction of the correction applied at each iteration.

4) Rescale globally: All contact counts in the dataset are multiplied by $\mathbf{f}$, ensuring that the adjustment is global and not distance-specific. 

5) After each scaling update, we re-fit the distance–decay curve of the rescaled dataset using the power-law form $$\mathbf c(s) = 10^{\,b + m \cdot \log_{10}(s)}$$ where $b$ is the intercept and $m$ is the slope obtained from linear regression in log–log space. The scaling exponent $\alpha = -m$ is preserved, so we adjust contact counts globally without artificially changing the decay rate of the dataset.

6) Measure the error: We then compare the fitted curve of the scaled dataset against the reference by computing a weighted mean squared error (MSE):

$$\text{WMSE} \;=\; \frac{\sum_{s} w(s) \, \bigl(\log c_{\text{scaled}}(s) - \log c_{\text{ref}}(s)\bigr)^2}{\sum_{s} w(s)}$$

where the weights $\mathbf{w(s) \propto s^{-1.08}}$ emphasize short genomic distances corrections. 
We use a weighted version of the error because short-range contacts dominate Hi-C datasets and are biologically more informative. Without weighting, imorivements at very large distances could overshadow meaningful corrections at shorter distances. By prioritizing these short distances, DifFracTion ensures that the normalization aligns most strongly with regions where biological signal is rich. 