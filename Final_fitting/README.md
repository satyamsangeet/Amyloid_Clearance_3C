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
\[
$\text{WRMSE} = \sqrt{\frac{\sum_{i=1}^{n} w_i \cdot (y_i - \hat{y}_i)^2}{\sum_{i=1}^{n} w_i}}$
\]

where:
- \$w_i\$ is the weight assigned to the \$i\$-th observation.
- \$y_i\$ are the true values.
- \$\hat{y}_i\$ are the predicted values.

- We have used NRMSE in the case of model fit with Huang et al., 2012 data, since the experimnetal data point doesn't have any error bar associated with it. In this case, all data points were treated equally. The optimization algorithm aimed to minimize the overall error across all points, without giving preference to any particular data points based on their uncertainty.

- We have used WRMSE in the case of model fit with Liu et al., 2022 data, since the experimnetal data point have error bars associated with it. In this case, we gave more weight to data points with lower uncertainty (low SD) and less weight to those with higher uncertainty (high SD). This means the optimization algorithm prioritized fitting the model to the more reliable data points, resulting in parameter values that best fit those specific points.
