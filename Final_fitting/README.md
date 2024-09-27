## Error Metrics

### Normalized RMSE (NRMSE)
Normalized RMSE (NRMSE) accounts for the scale of the target data, making it useful for comparing models across different datasets.

**Formula:**
\[
$\text{NRMSE} = \frac{\text{RMSE}}{\text{range}(y)}$
\]
or
\[
$\text{NRMSE} = \frac{\text{RMSE}}{\mu_y}$
\]

where:
- \$\text{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^{n} (y_i - \hat{y}_i)^2}$
- \$y_i\$ are the true values.
- \$\hat{y}_i\$ are the predicted values.
- \$\text{range}(y)\$ is the difference between the maximum and minimum values of \(y\).
- \$\mu_y\$ is the mean of the true values.

### Weighted RMSE (WRMSE)
Weighted RMSE (WRMSE) allows for differential weighting of errors, making it useful when some observations are more important than others.

**Formula:**
$`WRMSE = \sqrt{\frac{\sum_{i=1}^{n}w_{i}.\left( y_{i} - \hat{y}_i\right)^2}{\sum_{i=1}^{n}w_{i}}}`$

where:
- \$w_i\$ is the weight assigned to the \$i\$-th observation.
- \$y_i\$ are the true values.
- \$\hat{y}_i\$ are the predicted values.

**Note**

We have used NRMSE in the case of model fit with Huang et al., 2012 data, since the experimnetal data point doesn't have any error bar associated with it. In this case, all data points were treated equally. The optimization algorithm aimed to minimize the overall error across all points, without giving preference to any particular data points based on their uncertainty.

We have used WRMSE in the case of model fit with Liu et al., 2022 data, since the experimnetal data point have error bars associated with it. In this case, we gave more weight to data points with lower uncertainty (low SD) and less weight to those with higher uncertainty (high SD). This means the optimization algorithm prioritized fitting the model to the more reliable data points, resulting in parameter values that best fit those specific points.

**Optimised Model Parameters**

## Three Compartment Model (Brain, CSF, Plasma)

| **Parameter** | **Liu Fit Values** | **Huang Fit Values** | **Combined Fit Values** |
|---------------|-------------------|----------------------|----------------------|
| a             | 1.000010           | 2.945697       | 1.49 |
| b             | 1.603175           | 4.389072       | 1.49 |
| c             | 1.532420           | 4.006562        | 1.47 |
| A_wake (pg/ml/hr)        | 22.056825          | 21.074357       | 23.5 |
| A_sleep (pg/ml/hr)       | 0.8 * A_wake       | 0.8 * A_wake       | 0.8 * A_wake |
| a12_wake (/hr)      | 0.999998           | 0.116481        | 0.99 |
| a12_sleep (/hr)     | 2.5 * a12_wake     | 2.5 * a12_wake       | 2.5 * a12_wake |
| a13_wake (/hr)      | 0.100000           | 0.010355        | 0.09 |
| a13_sleep (/hr)     | a * a13_wake       | a * a13_wake       | a * a13_wake |
| a23_sleep (/hr)     | 0.06601            | 0.06601       | 0.06601 |
| a23_wake (/hr)      | a23_sleep / b      | a23_sleep / b       | a23_sleep / b |
| k_wake (/hr)        | 0.231049           | 0.231049       | 0.277258 |
| k_sleep (/hr)       | c * k_wake         | c * k_wake       | c * k_wake |

**Combined Fit** represents the model fitting on the normal sleep-wake data of Lucey2018 (CSF Data), Blattner2020 (CSF Data), Liu2022 (CSF + Plasma Data). WRMSE was used to fit the model with the experimnetal data as all the experimnetal data had stadard error associated with it
