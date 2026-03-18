# HPRU_HZQoL

Meta analysis for health-related quality of life loss

The repo performed a quantitative synthesis and statistical analysis for the HRQoL of herpes zoster patients .

## Project structure

### Modelling

-   **models/**: statistical models
-   **1\_**: exploring and testing models
-   **2\_**: model fitting
-   **posteriors**: simulated posterior parameters

### Simulation

-   **3_1_calculate_goodness_of_fit.R**: calculating and visualising the goodness of fit
-   **3_2_simulate_qol.R**: simulating quality of life
-   **3_3_simulate_shortfall.R**: simulating shortfall and quality-adjusted life expectancy

### Outputs

Fig a shows the distribution of HRQOl of HZ patients by time

We fristly test the assumptions for the proportional hypothesis. However, in this stage, we havn't estimates the event time for each data points. We used the data whose endpoints met and the proximate the event time by taking the mid point between the last noo-healthy measurement and the first healthy measurement for each patients.

Subjects were grouped into 10 age groups. We plotted the survival curve accordingly, finding the survival curves by age groups are not crossover to each other

With the proportional hazard hypothesis, we derived the respective likelihood functions with the time-to event distributions of WeiBull, expoenential, log-normal distributions. The prior distributions of the hyperparameters were set to align with the range of recovery time from the literatures of xx and xx days after rash onset. The posterior model were performed on systematicallu sampled 25% training samples, ensureing the samples from each study were evenly selected.

In summary, we observed.

1.  The age pattern of the time-to recovery can be captured through the proportional harzard model.
2.  The x distribution performed the best in terms the mean squrared errors but after considering thecompelexity and the overall goodness of fit, we decided to apply the model with eponential distributions on the time-to recovery data

## Quality of life after HZ rash onset

The temproal patterns were built on the (1) the proportional each states, (2) the values ofe HRQoL, and (3) both. For reducing the computation burden and there is clear boundary between clusters, Kmeans clustering were performed and data were labelled before exploring the age-time patterns.

From clinical perspectives, HZ were group by PHN and non-PHN according to a cut-off of 90 days after rash onset. PHN and can also be diagnosed if a patient surffering from extremely pain.

However, in our analysis were found that there is no specific level shift before and after 90 days. Instead, were not present differently according to fe

Patient-wise we also found that only x\$ of patients stayed in the same cluster through their disease course while the others frequently jump between severe-discomfort and mild discomfort states.

In the future it will be worth to do and study witch repeat measurement of HRQoL with equal time steps to reason the transitions and associated the transition to the treatment used and other covariates .

Population norm is defined as the baseline health status of a population. Because not all people in

However, the population norm estimates were usually overestimates since the institutionalised population and those who cannot repond the health survey questionaire were excluded.

The studies with the baseline health status surveyed were baseline the methods "then test". A retrospective survey of their health status after disease onset. The methods can detect the existance of response shift but is not able to specify the source of response shift.

## [License](LICENSE)
