# Connectome-based modeling of WISC full-scale IQ
Connectome-based predictive modeling (CPM)[^1][^2] is a validated supervised machine-learning approach. It is widely used to predict individual behavioral phenotypes, and social or clinical outcomes such as fluid intelligence[^2], sustained attention[^3][^4], personality traits[^5][^6], brain aging[^7], and abstinence[^8][^9].
CPM takes functional connectivity and one-dimensional behavioral data of each participant as input. Behavior-associated edges are selected based on the strength of their correlation (both positive and negative) with the behavior of interest by an arbitrary threshold[^10]. These edges are then summarized at the participant level, which results in a summary statistic, ‘network strength’, for each participant. The network strength and behavior data are used to fit a simple linear model where network strength is the independent variable and behavior score is the dependent variable. Thus, the network associated with the behavior of interest is built and can be used to predict other individuals' behavior from their functional connectivity data. Evaluation of the model is based on its predictive ability and statistical significance.

## CPM Python Adaptation
- Based on [CPM Tutorial](https://github.com/esfinn/cpm_tutorial).
- Cross-validation
    - Added seed for reproduction.
- Feature selecion
    - Added p-vlaue threshold.
    - Added partial correlation and robust regression.
    - Pandas apply was used in `cpm_pycpmFunctions.py`.
    - For-loop was used in `cpm_pycpmFunctionsLoop.py`.


## Key Points
1. Data normality.
2. Correlation between behavior and functional connectivity.
3. *k* in cross-validation.
4. Correlation method in feature seleciton.
5. Edges summarization/averaging.
5. Permutation testing for statistical significance.




[^1]: Finn, E. S., Shen, X., Scheinost, D., Rosenberg, M. D., Huang, J., Chun, M. M., Papademetris, X., & Constable, R. T. (2015). Functional connectome fingerprinting: Identifying individuals using patterns of brain connectivity. *Nature Neuroscience*, *18*(11), Article 11. [Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity | Nature Neuroscience](https://doi.org/10.1038/nn.4135)
[^2]: Shen, X., Finn, E. S., Scheinost, D., Rosenberg, M. D., Chun, M. M., Papademetris, X., & Constable, R. T. (2017). Using connectome-based predictive modeling to predict individual behavior from brain connectivity. _Nature Protocols_, _12_(3), Article 3. [https://doi.org/10.1038/nprot.2016.178](https://doi.org/10.1038/nprot.2016.178)
[^3]: Rosenberg, M. D., Scheinost, D., Greene, A. S., Avery, E. W., Kwon, Y. H., Finn, E. S., Ramani, R., Qiu, M., Constable, R. T., & Chun, M. M. (2020). Functional connectivity predicts changes in attention observed across minutes, days, and months. *Proceedings of the National Academy of Sciences*, *117*(7), 3797–3807. https://doi.org/10.1073/pnas.1912226117 
[^4]: Yoo, K., Rosenberg, M. D., Hsu, W.-T., Zhang, S., Li, C.-S. R., Scheinost, D., Constable, R. T., & Chun, M. M. (2018). Connectome-based predictive modeling of attention: Comparing different functional connectivity features and prediction methods across datasets. *NeuroImage*, *167*, 11–22. https://doi.org/10.1016/j.neuroimage.2017.11.010
[^5]: Cai, H., Zhu, J., & Yu, Y. (2020). Robust prediction of individual personality from brain functional connectome. *Social Cognitive and Affective Neuroscience*, *15*(3), 359–369. https://doi.org/10.1093/scan/nsaa044
[^6]: Dubois, J., Galdi, P., Han, Y., Paul, L. K., & Adolphs, R. (2018). Resting-State Functional Brain Connectivity Best Predicts the Personality Dimension of Openness to Experience. *Personality Neuroscience*, *1*, e6. https://doi.org/10.1017/pen.2018.8
[^7]: Kim, E., Kim, S., Kim, Y., Cha, H., Lee, H. J., Lee, T., & Chang, Y. (2022). Connectome-based predictive models using resting-state fMRI for studying brain aging. *Experimental Brain Research*. https://doi.org/10.1007/s00221-022-06430-7
[^8]: Lichenstein, S. D., Scheinost, D., Potenza, M. N., Carroll, K. M., & Yip, S. W. (2021). Dissociable neural substrates of opioid and cocaine use identified via connectome-based modelling. *Molecular Psychiatry*, *26*(8), 4383–4393. https://doi.org/10.1038/s41380-019-0586-y
[^9]: Yip, S. W., Scheinost, D., Potenza, M. N., & Carroll, K. M. (2019). Connectome-Based Prediction of Cocaine Abstinence. *American Journal of Psychiatry*, *176*(2), 156–164. https://doi.org/10.1176/appi.ajp.2018.17101147
[^10]: Greene, A. S., Gao, S., Scheinost, D., & Constable, R. T. (2018). Task-induced brain state manipulation improves prediction of individual traits. _Nature Communications_, _9_(1), 2807. [https://doi.org/10.1038/s41467-018-04920-3](https://doi.org/10.1038/s41467-018-04920-3)   