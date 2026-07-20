import numpy as np
import matplotlib.pyplot as plt

# params
r_bc = 0.038
r_bp = 0.014
r_cp = 0.00537
sigma_bc = 1.131
sigma_bp = 1.768
sigma_cp = 6.100
A = 16.203
sigma_A = 0.772
r_p = 0.427
sigma_p = 4.253

# params to analyze
parameters_config = {
    'sigma_p': {
        'values': np.linspace(0.5, 7, 30),
        'default': sigma_p,
        'vertical_line': sigma_p,
        'label': 'σ_p'
    },
    'r_p': {
        'values': np.linspace(0.3, 0.6, 30),
        'default': r_p,
        'vertical_line': r_p, 
        'label': 'r_p'
    },
    'A': {
        'values': np.linspace(0, 80, 30),
        'default': A,
        'vertical_line': A,
        'label': 'A'
    },
    'sigma_A': {
        'values': np.linspace(0, 0.99, 30),
        'default': sigma_A,
        'vertical_line': sigma_A,
        'label': 'σ_A'
    },
    'r_bc': {
        'values': np.linspace(0, 0.5, 30),
        'default': r_bc,
        'vertical_line': r_bc,
        'label': 'r_bc'
    },
    'sigma_bc': {
        'values': np.linspace(0.1, 7, 30),
        'default': sigma_bc,
        'vertical_line': sigma_bc,
        'label': 'σ_bc'
    },
    'r_bp': {
        'values': np.linspace(0, 0.5, 30),
        'default': r_bp,
        'vertical_line': r_bp,
        'label': 'r_bp'
    },
    'sigma_bp': {
        'values': np.linspace(0.1, 7, 30),
        'default': sigma_bp,
        'vertical_line': sigma_bp,
        'label': 'σ_bp'
    },
    'r_cp': {
        'values': np.linspace(0, 0.06, 30),
        'default': r_cp,
        'vertical_line': r_cp,
        'label': 'r_cp'
    },
    'sigma_cp': {
        'values': np.linspace(0.3, 7, 30),
        'default': sigma_cp,
        'vertical_line': sigma_cp,
        'label': 'σ_cp'
    }
}

def calculate_steady_states(param_name, param_values):
    brain_wake = np.zeros_like(param_values)
    brain_sleep = np.zeros_like(param_values)
    csf_wake = np.zeros_like(param_values)
    csf_sleep = np.zeros_like(param_values)
    plasma_wake = np.zeros_like(param_values)
    plasma_sleep = np.zeros_like(param_values)
    
    for i, param_val in enumerate(param_values):
        current_r_bc = r_bc if param_name != 'r_bc' else param_val
        current_r_bp = r_bp if param_name != 'r_bp' else param_val
        current_r_cp = r_cp if param_name != 'r_cp' else param_val
        current_sigma_bc = sigma_bc if param_name != 'sigma_bc' else param_val
        current_sigma_bp = sigma_bp if param_name != 'sigma_bp' else param_val
        current_sigma_cp = sigma_cp if param_name != 'sigma_cp' else param_val
        current_A = A if param_name != 'A' else param_val
        current_sigma_A = sigma_A if param_name != 'sigma_A' else param_val
        current_r_p = r_p if param_name != 'r_p' else param_val
        current_sigma_p = parameters_config['sigma_p']['default'] if param_name != 'sigma_p' else param_val
        
        # Calculate steady states
        brain_wake[i] = current_A / (current_r_bc + current_r_bp)
        brain_sleep[i] = (current_sigma_A * current_A) / (current_sigma_bc * current_r_bc + current_sigma_bp * current_r_bp)
        
        csf_wake[i] = (current_r_bc * current_A) / (current_r_cp * (current_r_bc + current_r_bp))
        csf_sleep[i] = (current_sigma_bc * current_r_bc * current_sigma_A * current_A) / (current_sigma_cp * current_r_cp * (current_sigma_bc * current_r_bc + current_sigma_bp * current_r_bp))
        
        plasma_wake[i] = current_A / current_r_p
        plasma_sleep[i] = (current_sigma_A * current_A) / (current_sigma_p * current_r_p)
    
    return brain_wake, brain_sleep, csf_wake, csf_sleep, plasma_wake, plasma_sleep

def create_plots(param_name, param_values, brain_wake, brain_sleep, csf_wake, csf_sleep, plasma_wake, plasma_sleep, param_config):
    param_label = param_config['label']
    vertical_line = param_config['vertical_line']
    
    # Brain
    plt.figure(figsize=(6,6))
    plt.plot(param_values, brain_wake, label='Wake', marker="o", markersize=8, linewidth=3.0)
    plt.plot(param_values, brain_sleep, label='Sleep', marker="o", markersize=8, linewidth=3.0)
    plt.axvline(x=vertical_line, linestyle="--", color="k", linewidth=2, alpha=0.7)
    plt.grid(alpha=0.1, color="gray", linewidth=2.5, linestyle="--")
    plt.xlabel(param_label, fontsize=15)
    plt.ylabel('Brain Concentration', fontsize=15)
    plt.title(f'Brain vs {param_label}', fontsize=16)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.savefig(f"brain_with_{param_name}.png", dpi=300, bbox_inches="tight")
    plt.show()
    
    # CSF
    plt.figure(figsize=(6,6))
    plt.plot(param_values, csf_wake, label='Wake', marker="o", markersize=8, linewidth=3.0)
    plt.plot(param_values, csf_sleep, label='Sleep', marker="o", markersize=8, linewidth=3.0)
    plt.axhspan(500, 1000, alpha=0.2, color='green')
    plt.axvline(x=vertical_line, linestyle="--", color="k", linewidth=2, alpha=0.7)
    plt.grid(alpha=0.1, color="gray", linewidth=2.5, linestyle="--")
    plt.xlabel(param_label, fontsize=15)
    plt.ylabel('CSF Concentration', fontsize=15)
    plt.title(f'CSF vs {param_label}', fontsize=16)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.savefig(f"csf_with_{param_name}.png", dpi=300, bbox_inches="tight")
    plt.show()
    
    # Plasma
    plt.figure(figsize=(6,6))
    plt.plot(param_values, plasma_wake, label='Wake', marker="o", markersize=8, linewidth=3.0)
    plt.plot(param_values, plasma_sleep, label='Sleep', marker="o", markersize=8, linewidth=3.0)
    plt.axhspan(15, 20, alpha=0.2, color='green')
    plt.axvline(x=vertical_line, linestyle="--", color="k", linewidth=2, alpha=0.7)
    plt.grid(alpha=0.1, color="gray", linewidth=2.5, linestyle="--")
    plt.xlabel(param_label, fontsize=15)
    plt.ylabel('Plasma Concentration', fontsize=15)
    plt.title(f'Plasma vs {param_label}', fontsize=16)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.savefig(f"plasma_with_{param_name}.png", dpi=300, bbox_inches="tight")
    plt.show()

# Main
def run_analysis():
    for param_name, param_config in parameters_config.items():
        print(f"\nAnalyzing parameter: {param_name}")
        param_values = param_config['values']

        brain_wake, brain_sleep, csf_wake, csf_sleep, plasma_wake, plasma_sleep = calculate_steady_states(param_name, param_values)
        create_plots(param_name, param_values, brain_wake, brain_sleep, csf_wake, csf_sleep, plasma_wake, plasma_sleep, param_config)
        print(f"Plots saved for {param_name}")

# Run
if __name__ == "__main__":
    run_analysis()
    print("\nAnalysis complete! All plots have been generated and saved.")
