import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import pandas as pd
from tqdm import tqdm

# Read experimental data from CSV files
data_exp1 = np.loadtxt('/content/huang_csf42.csv', delimiter=',', skiprows=1)
data_exp2 = np.loadtxt('/content/huang_plasma42.csv', delimiter=',', skiprows=1)
time_exp = data_exp1[:, 0]
exp_csf = data_exp1[:, 1]
exp_plasma = data_exp2[:, 1]

def model(t, y, a21, a12_wake, a12_sleep, A_wake, A_sleep, k):
    sw_cycle = (t % 24 >= 8) & (t % 24 < 24)
    dydt_n = np.zeros(2)
    dydt_n[0] = A_wake * sw_cycle + A_sleep * (1 - sw_cycle) + a21 * y[1] - (a12_wake * sw_cycle + a12_sleep * (1 - sw_cycle)) * y[0]
    dydt_n[1] = (a12_wake * sw_cycle + a12_sleep * (1 - sw_cycle)) * y[0] - (a21 + k) * y[1]
    return dydt_n

def objective_function(params, exp_csf, exp_plasma):
    a21, a12_wake, a12_sleep, A_wake, A_sleep, k = params
    sol = solve_ivp(lambda t, y: model(t, y, a21, a12_wake, a12_sleep, A_wake, A_sleep, k), [0, 24*10], [exp_csf[0], exp_plasma[0]], max_step=0.01)

    csf_last_36hours_data = sol.y[0, 20000:23601]
    plasma_last_36hours_data = sol.y[1, 20000:23601]

    time_indices = range(1, 37)
    csf_last_36hours = csf_last_36hours_data[time_indices]
    plasma_last_36hours = plasma_last_36hours_data[time_indices]

    avg_C = np.mean(csf_last_36hours[:12])
    avg_B = np.mean(plasma_last_36hours[:12])

    normalise_C = (csf_last_36hours - avg_C) / avg_C
    normalise_B = (plasma_last_36hours - avg_B) / avg_B

    error_C = np.sqrt(np.mean((normalise_C - exp_csf) ** 2))
    error_B = np.sqrt(np.mean((normalise_B - exp_plasma) **2))

    total_error = error_C + error_B
    return total_error

# Define initial guess for parameters
initial_guess = [0.1, 0.1, 0.3, 30, 1, 2]

# Define callback function for tqdm
def callback(x):
    pbar.update()

# Store loss values
loss_values = []

# Minimize the objective function to find optimal parameters
with tqdm(total=50) as pbar:
    def callback_loss(x):
        loss_values.append(objective_function(x, exp_csf, exp_plasma))
        callback(x)

    result = minimize(objective_function, initial_guess, args=(exp_csf, exp_plasma), bounds=((0, 1), (0, 1), (0, 1), (20, 50), (0, 10), (0, 10)), callback=callback_loss)

# Get optimized parameters
a21_opt, a12_wake_opt, a12_sleep_opt, A_wake_opt, A_sleep_opt, k_opt = result.x

# Print optimized parameters
print('Optimized Parameters:')
print('A_wake:', A_wake_opt)
print('A_sleep:', A_sleep_opt)
print('a12_wake:', a12_wake_opt)
print('a12_sleep:', a12_sleep_opt)
print('a21:', a21_opt)
print('k:', k_opt)

# Solve the model with optimized parameters at each hour and save the results
results = []
for hour in tqdm(range(241)):
    sol = solve_ivp(lambda t, y: model(t, y, a21_opt, a12_wake_opt, a12_sleep_opt, A_wake_opt, A_sleep_opt, k_opt), [0, hour], [exp_csf[0], exp_plasma[0]], max_step=0.01)
    results.append([hour, sol.y[0][-1], sol.y[1][-1]])

# Save results to CSV
df = pd.DataFrame(results, columns=['Time', 'CSF', 'Plasma'])
df.to_csv('/content/model_v20_output_py.csv', index=False)

# Save loss values to CSV
loss_df = pd.DataFrame(loss_values, columns=['Loss'])
loss_df.to_csv('/content/loss_values.csv', index=False)
