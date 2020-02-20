

# Sensitivity of optimal choice of vaccination strategy to changes in social contact patterns
## Abstract

*Background* The Bayesian framework allows one to use prior information when calibrating transmission models of vaccination policy. Previous transmission models have not assessed the value of prior information of the contact patterns. However, many vaccine policy decisions are based on quantifying the value of indirect effects from vaccination, such as LAIV. Many countries do not have empirical studies describing their contact patterns, and often those of geographical neighbours are used.<br/>
*Methods* We fit and parameterized a transmission model of influenza in Great Britain using survey data from 10 other countries to construct counterfactual contact patterns. Using the counterfactual contact patterns we conducted cost-effectiveness analysis of three pediatric vaccination policies and compared the results to those derived using the country-specific Great Britain contact structure.<br/>
*Results* Country-specific contact patterns add value in calibrating transmission models. Specifically, the prior information on contact patterns do not converge on the same posterior, such that one’s posterior belief about the relative importance of each age group in transmission changes based on the prior information. We show the relative impact of age groups on transmission that come from the choice of contact patterns used as priors decreases the probability a strategy is considered cost-effective by propagating uncertainty through the model.<br/>
*Conclusion* Even using substitute patterns as prior information can alter vaccine decision-making. Therefore, we would advocate that in the absence of country-specific contact patterns, structural uncertainty analysis should be undertaken when quantifying the effect of strategies that involve herd immunity.

### Data Sources
We obtained contact data from multiple sources most notably the POLYMOD study (Mossong et al., 2008), but also from similar studies in Zimbabwe, Peru, and France using the [Zenodo database](https://zenodo.org/communities/social_contact_data?page=1&size=20). Raw survey data  was downloaded and weighted by contacts for weekday vs weekend following the method described by Baguelin et al. (2015) using the code in [INPUT_polymod_pull](https://github.com/tajwenzel/IntlSens/blob/master/INPUT_polymod_pull.R). 

#### Other Data

Other data which was needed for these analyses are included in the repository:
* Vaccine uptake rate functions for final coverage estimated at 30%, 55%, and 70% starting on September 1st and continuing till December 12th in [cov.function30.Rdata](https://github.com/tajwenzel/IntlSens/blob/master/cov.function30.RData), [cov.function55.Rdata](https://github.com/tajwenzel/IntlSens/blob/master/cov.function55.RData), [cov.function70.Rdata](https://github.com/tajwenzel/IntlSens/blob/master/cov.function70.RData), respectively.
* New LAIV and TIIV vaccine costs by delivery method laid out in Supplemental Section, Table 5.
* Prior distributions for Markov chain Monte Carlo laid out in Supplemental Section Table 3.
* Incidence of Influenza-like-illnesses from 1995-2009 in the UK obtained from the 'fluEvidenceSynthesis' R package. 

