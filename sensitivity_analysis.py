import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd

class Model:
    def __init__(self):
        self.params = {
            'r_bc': 0.038,
            'r_bp': 0.014,
            'r_cp': 0.00537,
            'sigma_bc': 1.131,
            'sigma_bp': 1.768,
            'sigma_cp': 6.100,
            'sigma_p': 4.253,
            'A': 16.203,
            'sigma_A': 0.772,
            'r_p': 0.427
        }
    
    def model_equations(self, t, y, params=None):
        if params is None:
            params = self.params

        r_bc = params['r_bc']
        r_bp = params['r_bp']
        r_cp = params['r_cp']
        sigma_bc = params['sigma_bc']
        sigma_bp = params['sigma_bp']
        sigma_cp = params['sigma_cp']
        sigma_p = params['sigma_p']
        A = params['A']
        sigma_A = params['sigma_A']
        r_p = params['r_p']
        
        # Switch
        sw_cycle = 1 if (t % 24 >= 8 and t % 24 < 24) else 0
        dydt = np.zeros(3)
        
        # Brain
        dydt[0] = (A * sw_cycle + sigma_A * A * (1 - sw_cycle) - 
                   (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + 
                    r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y[0])
        
        # CSF
        dydt[1] = ((r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y[0] - 
                   (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y[1])
        
        # Plasma
        dydt[2] = ((r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y[0] + 
                   (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y[1] - 
                   (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y[2])
        
        return dydt
    
    def simulate(self, t_span=(0, 24), y0=None, params=None):
        if y0 is None:
            y0 = [0, 0, 0]  #IC
        
        if params is None:
            params = self.params

        sol = solve_ivp(
            lambda t, y: self.model_equations(t, y, params),
            t_span, y0, 
            dense_output=True,
            rtol=1e-8, atol=1e-10
        )
        
        return sol

def perform_sensitivity_analysis(model, perturbation=0.01, stabilization_time=2400, 
                               analysis_start=2330, analysis_end=2353, y0=None, output_dir='./', 
                               integration_method='trapz'):
    
    if y0 is None:
        y0 = [0, 0, 0]
    
    # Step 1
    print("Step 1: Running simulation to reach steady state...")
    print(f"Simulating from t=0 to t={stabilization_time} hours...")
    baseline_sol = model.simulate((0, stabilization_time), y0)
    
    # Extract data
    print(f"Step 2: Extracting data from analysis window: t={analysis_start} to t={analysis_end} hours")
    t_analysis = np.linspace(analysis_start, analysis_end, 1000)
    baseline_y_analysis = baseline_sol.sol(t_analysis)
    
    # system stabaliuse
    t_first_cycle = np.linspace(0, 24, 1000)
    t_last_cycle = np.linspace(stabilization_time-24, stabilization_time, 1000)
    baseline_y_first = baseline_sol.sol(t_first_cycle)
    baseline_y_last = baseline_sol.sol(t_last_cycle)
    
    print("Stabilization check:")
    for i, comp in enumerate(['Brain', 'CSF', 'Plasma']):
        first_mean = np.mean(baseline_y_first[i, :])
        last_mean = np.mean(baseline_y_last[i, :])
        change_percent = abs((last_mean - first_mean) / first_mean * 100) if first_mean != 0 else 0
        print(f"  {comp}: First cycle mean = {first_mean:.6f}, Last cycle mean = {last_mean:.6f}")
        print(f"         Change: {change_percent:.4f}% (should be <1% for good stabilization)")
    print()
    
    # Step 3
    print("Step 3: Integrating compartment concentrations...")
    if integration_method == 'trapz':
        baseline_outputs = np.trapz(baseline_y_analysis, t_analysis, axis=1)
    else:
        baseline_outputs = np.sum(baseline_y_analysis, axis=1) * (t_analysis[1] - t_analysis[0])
    
    print(f"Baseline integrated concentrations (analysis window):")
    for i, comp in enumerate(['Brain', 'CSF', 'Plasma']):
        print(f"  {comp}: {baseline_outputs[i]:.8f}")
    print()
    
    param_names = list(model.params.keys())
    compartment_names = ['Brain', 'CSF', 'Plasma']
    n_params = len(param_names)
    n_compartments = len(compartment_names)
    
    sensitivity_matrix = np.zeros((n_compartments, n_params))
    perturbed_outputs = {}
    
    # Step 4: Perturb 1%
    print("Step 4: Perturbing each parameter independently by 1%...")
    print()
    
    for j, param_name in enumerate(param_names):
        print(f"Analyzing parameter {j+1}/{n_params}: {param_name} = {model.params[param_name]}")
        
        perturbed_params = model.params.copy()
        original_value = model.params[param_name]
        perturbed_value = original_value * (1 + perturbation)
        perturbed_params[param_name] = perturbed_value
        
        print(f"  Original: {original_value:.6f}")
        print(f"  Perturbed: {perturbed_value:.6f}")

        perturbed_sol = model.simulate((0, stabilization_time), y0, perturbed_params)
        perturbed_y_analysis = perturbed_sol.sol(t_analysis)

        if integration_method == 'trapz':
            perturbed_integrated = np.trapz(perturbed_y_analysis, t_analysis, axis=1)
        else:
            perturbed_integrated = np.sum(perturbed_y_analysis, axis=1) * (t_analysis[1] - t_analysis[0])
        perturbed_outputs[param_name] = perturbed_integrated
        
        # Step 5: Calculate normalized sensitivity
        for i in range(n_compartments):
            delta_y_i = perturbed_integrated[i] - baseline_outputs[i]
            delta_p_j = perturbed_value - original_value
            
            if baseline_outputs[i] != 0 and original_value != 0:
                sensitivity_matrix[i, j] = (delta_y_i / baseline_outputs[i]) / (delta_p_j / original_value)
                print(f"  {compartment_names[i]}: S = {sensitivity_matrix[i, j]:.6f}")
            else:
                sensitivity_matrix[i, j] = 0
                print(f"  {compartment_names[i]}: S = 0 (zero baseline)")
        print()
    
    detailed_results = {
        'baseline_outputs': baseline_outputs,
        'perturbed_outputs': perturbed_outputs,
        'sensitivity_matrix': sensitivity_matrix,
        'param_names': param_names,
        'compartment_names': compartment_names,
        'parameters': model.params,
        'perturbation': perturbation,
        'stabilization_time': stabilization_time,
        'analysis_window': (analysis_start, analysis_end)
    }
  
    save_results_to_text(detailed_results, output_dir)
    return sensitivity_matrix, baseline_outputs, detailed_results

def save_results_to_text(results, output_dir='./'):
    with open(f'{output_dir}sensitivity_analysis_report.txt', 'w') as f:
        f.write("LOCAL SENSITIVITY ANALYSIS REPORT\n")
        f.write("="*50 + "\n\n")
        
        f.write("METHODOLOGY:\n")
        f.write(f"- Stabilization time: {results['stabilization_time']} hours ({results['stabilization_time']/24:.1f} days)\n")
        f.write(f"- Analysis window: {results['analysis_window'][0]} to {results['analysis_window'][1]} hours\n")
        f.write(f"- Parameter perturbation: {results['perturbation']*100}%\n")
        f.write("- Each parameter increased by 1% while others held constant\n")
        f.write("- Long simulation to reach steady state, then analysis on final 24-hour cycle\n")
        f.write("- Concentrations integrated over analysis window\n")
        f.write("- Normalized sensitivity: S_ij = (Δy_i/y_i) / (Δp_j/p_j)\n\n")
        
        f.write("BASELINE PARAMETERS:\n")
        for param, value in results['parameters'].items():
            f.write(f"{param:10s} = {value:12.6f}\n")
        f.write("\n")
        
        f.write("BASELINE INTEGRATED CONCENTRATIONS (final 24-hour cycle):\n")
        for i, compartment in enumerate(results['compartment_names']):
            f.write(f"{compartment:10s} = {results['baseline_outputs'][i]:15.6f}\n")
        f.write("\n")
        
        f.write("NORMALIZED SENSITIVITY COEFFICIENTS:\n")
        f.write("Compartment".ljust(12))
        for param in results['param_names']:
            f.write(param.ljust(12))
        f.write("\n")
        f.write("-" * (12 + 12*len(results['param_names'])) + "\n")
        
        for i, compartment in enumerate(results['compartment_names']):
            f.write(compartment.ljust(12))
            for j in range(len(results['param_names'])):
                f.write(f"{results['sensitivity_matrix'][i,j]:11.6f} ")
            f.write("\n")
        f.write("\n")
        
        f.write("MOST SENSITIVE PARAMETERS (by absolute value):\n")
        for i, compartment in enumerate(results['compartment_names']):
            most_sensitive_idx = np.argmax(np.abs(results['sensitivity_matrix'][i, :]))
            most_sensitive_param = results['param_names'][most_sensitive_idx]
            sensitivity_value = results['sensitivity_matrix'][i, most_sensitive_idx]
            f.write(f"{compartment:10s}: {most_sensitive_param:10s} (S = {sensitivity_value:8.4f})\n")
    
    # 2. Sensitivity matrix
    with open(f'{output_dir}sensitivity_matrix.txt', 'w') as f:
        f.write("# Normalized Sensitivity Coefficients Matrix\n")
        f.write("# Rows: Compartments, Columns: Parameters\n")
        f.write("# S_ij = (delta_y_i/y_i) / (delta_p_j/p_j)\n")
        f.write("Compartment,")
        f.write(",".join(results['param_names']) + "\n")
        
        for i, compartment in enumerate(results['compartment_names']):
            f.write(f"{compartment},")
            f.write(",".join([f"{results['sensitivity_matrix'][i,j]:.8f}" for j in range(len(results['param_names']))]) + "\n")
    
    # 3. Detailed parameter information
    with open(f'{output_dir}parameter_details.txt', 'w') as f:
        f.write("PARAMETER DETAILS AND PERTURBATION ANALYSIS\n")
        f.write("="*50 + "\n\n")
        
        for j, param_name in enumerate(results['param_names']):
            f.write(f"PARAMETER: {param_name}\n")
            f.write(f"Baseline value: {results['parameters'][param_name]:.8f}\n")
            f.write(f"Perturbed value: {results['parameters'][param_name] * (1 + results['perturbation']):.8f}\n")
            f.write(f"Absolute change: {results['parameters'][param_name] * results['perturbation']:.8f}\n")
            f.write(f"Relative change: {results['perturbation']*100:.2f}%\n")
            f.write("\nSensitivity coefficients:\n")
            for i, compartment in enumerate(results['compartment_names']):
                f.write(f"  {compartment:10s}: {results['sensitivity_matrix'][i,j]:10.6f}\n")
            f.write("\n" + "-"*40 + "\n\n")
    
    # 4. Baseline and perturbed outputs
    with open(f'{output_dir}concentration_outputs.txt', 'w') as f:
        f.write("INTEGRATED CONCENTRATIONS (24-hour sum)\n")
        f.write("="*45 + "\n\n")
        
        f.write("BASELINE OUTPUTS:\n")
        for i, compartment in enumerate(results['compartment_names']):
            f.write(f"{compartment:10s}: {results['baseline_outputs'][i]:15.8f}\n")
        f.write("\n")
        
        f.write("PERTURBED OUTPUTS (for each parameter perturbation):\n")
        for param_name in results['param_names']:
            f.write(f"\nParameter: {param_name}\n")
            perturbed_vals = results['perturbed_outputs'][param_name]
            for i, compartment in enumerate(results['compartment_names']):
                change = perturbed_vals[i] - results['baseline_outputs'][i]
                percent_change = (change / results['baseline_outputs'][i]) * 100 if results['baseline_outputs'][i] != 0 else 0
                f.write(f"  {compartment:10s}: {perturbed_vals[i]:15.8f} (change: {change:+10.8f}, {percent_change:+7.3f}%)\n")
    
    print(f"✓ Summary report: {output_dir}sensitivity_analysis_report.txt")
    print(f"✓ Sensitivity matrix: {output_dir}sensitivity_matrix.txt")
    print(f"✓ Parameter details: {output_dir}parameter_details.txt")
    print(f"✓ Concentration outputs: {output_dir}concentration_outputs.txt")
    print()

def plot_sensitivity_results(sensitivity_matrix, param_names, compartment_names, save_plots=True):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Heatmap
    im = ax1.imshow(sensitivity_matrix, cmap='RdBu_r', aspect='auto')
    ax1.set_xticks(range(len(param_names)))
    ax1.set_xticklabels(param_names, rotation=45, ha='right')
    ax1.set_yticks(range(len(compartment_names)))
    ax1.set_yticklabels(compartment_names)
    ax1.set_title('Normalized Sensitivity Coefficients\nS_ij = (Δy_i/y_i) / (Δp_j/p_j)')

    plt.colorbar(im, ax=ax1, label='Sensitivity Coefficient')

    for i in range(len(compartment_names)):
        for j in range(len(param_names)):
            text = ax1.text(j, i, f'{sensitivity_matrix[i, j]:.3f}',
                           ha="center", va="center", color="black", fontsize=9)

    abs_sensitivity = np.abs(sensitivity_matrix)
    x_pos = np.arange(len(param_names))
    width = 0.25
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    for i, compartment in enumerate(compartment_names):
        ax2.bar(x_pos + i*width, abs_sensitivity[i, :], width, 
                label=compartment, alpha=0.8, color=colors[i])
    
    ax2.set_xlabel('Parameters')
    ax2.set_ylabel('|Sensitivity Coefficient|')
    ax2.set_title('Absolute Sensitivity by Parameter')
    ax2.set_xticks(x_pos + width)
    ax2.set_xticklabels(param_names, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plots:
        plt.savefig('sensitivity_analysis_plots.png', dpi=300, bbox_inches='tight')
        plt.savefig('sensitivity_analysis_plots.pdf', bbox_inches='tight')
        print("✓ Plots saved: sensitivity_analysis_plots.png and .pdf")
    
    plt.show()

# Main
if __name__ == "__main__":
    print("AMYLOID-β CONCENTRATION MODEL - LOCAL SENSITIVITY ANALYSIS")
    print("Following exact methodology from writeup")
    print()

    model = Model()
    print("Model parameters (globally optimized):")
    for param, value in model.params.items():
        print(f"  {param:10s} = {value:10.6f}")
    print()
    
    y0 = [0, 0, 0]  # IC

    print("AMYLOID-β CONCENTRATION MODEL - LOCAL SENSITIVITY ANALYSIS")
    print("Using long stabilization protocol as requested")
    print()
    
    sensitivity_matrix, baseline_outputs, detailed_results = perform_sensitivity_analysis(
        model, 
        perturbation=0.01, 
        stabilization_time=2400,  # 100 days 
        analysis_start=2336,      # Start of final 24-hour cycle
        analysis_end=2360,        # End of final 24-hour cycle  
        y0=y0
    )

    print("="*60)
    print("KEY RESULTS SUMMARY")
    print("="*60)
    
    print(f"\nBaseline integrated concentrations (final 24-hour cycle):")
    for i, compartment in enumerate(detailed_results['compartment_names']):
        print(f"  {compartment:10s}: {baseline_outputs[i]:15.8f}")
    
    print(f"\nMost sensitive parameters (highest |S_ij|):")
    for i, compartment in enumerate(detailed_results['compartment_names']):
        most_sensitive_idx = np.argmax(np.abs(sensitivity_matrix[i, :]))
        most_sensitive_param = detailed_results['param_names'][most_sensitive_idx]
        sensitivity_value = sensitivity_matrix[i, most_sensitive_idx]
        print(f"  {compartment:10s}: {most_sensitive_param:10s} (S = {sensitivity_value:8.4f})")
    
    plot_sensitivity_results(sensitivity_matrix, detailed_results['param_names'], 
                           detailed_results['compartment_names'])
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!")
    print("="*60)
    print("Long stabilization protocol used:")
    print(f"- Simulated for {detailed_results['stabilization_time']} hours ({detailed_results['stabilization_time']/24:.0f} days)")
    print(f"- Analysis performed on final 24-hour cycle (t={detailed_results['analysis_window'][0]} to {detailed_results['analysis_window'][1]} hours)")
    print("- This should provide much more stable and representative sensitivity coefficients")
    print("- All results saved to text files for documentation")
