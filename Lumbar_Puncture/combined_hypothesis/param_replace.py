import re
import os
from typing import Dict

def update_matlab_parameters(file_path: str, new_params: Dict[str, float], file_type: str) -> bool:
    """
    Updates parameter values in a MATLAB script file.
    
    Args:
        file_path: Path to the input MATLAB file
        new_params: Dictionary of parameter names and their new values
        file_type: Type of file (blattner, lucey, or liu) for logging purposes
        
    Returns:
        bool: True if successful, False otherwise
    """
    
    # Read the file
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return False
    except Exception as e:
        print(f"Error reading file '{file_path}': {e}")
        return False
    
    print(f"\nUpdating {file_type} parameters in: {file_path}")
    print("-" * 60)
    
    # Track if any parameters were actually updated
    updated_count = 0
    
    # Update each parameter
    for param_name, new_value in new_params.items():
        # Match decimals or scientific notation
        pattern = rf'({param_name}\s*=\s*)([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?);'
        
        matches = list(re.finditer(pattern, content))
        if matches:
            old_values = [m.group(2) for m in matches]
            content = re.sub(pattern, rf'\g<1>{new_value};', content)
            print(f"  {param_name}: {old_values[0]} -> {new_value}")
            updated_count += 1
    
    # Handle bounds updates
    bounds_rules = {
        'sigma_A': {'lb': new_params['sigma_A'], 'ub': 0.99},
        'sigma_bc': {'lb': 0, 'ub': new_params['sigma_bc']},
        'sigma_bp': {'lb': 0, 'ub': new_params['sigma_bp']},
        'sigma_cp': {'lb': 0, 'ub': new_params['sigma_cp']},
        'sigma_p': {'lb': 0, 'ub': new_params['sigma_p']},
        'r_bc': {'lb': new_params['r_bc'], 'ub': 5},
        'r_bp': {'lb': 0, 'ub': 0.25},
        'r_cp': {'lb': 0, 'ub': new_params['r_cp']},
        'r_p':  {'lb': 0, 'ub': 0.6},
    }
    
    def detect_optimized_parameters(content):
        param_pattern = r'([A-Za-z_]\w*)\s*=\s*params\((\d+)\)'
        matches = re.findall(param_pattern, content)
        
        if matches:
            sorted_matches = sorted(matches, key=lambda x: int(x[1]))
            seen, unique = set(), []
            for name, idx in sorted_matches:
                if idx not in seen:
                    unique.append((name, idx))
                    seen.add(idx)
            param_order = [m[0] for m in unique]
            return param_order
        else:
            return None
    
    def update_bounds_array(content, bound_type, bounds_rules, param_order):
        if not param_order:
            return content
        
        pattern = rf'({bound_type}\s*=\s*\[)([^\]]+)(\];)'
        matches = list(re.finditer(pattern, content))
        
        if not matches:
            return content
        
        for match in matches:
            bounds_str = match.group(2)
            bounds_values = [v.strip() for v in bounds_str.split(',')]
            
            if len(bounds_values) != len(param_order):
                continue
            
            new_values = []
            bounds_updated = False
            for i, name in enumerate(param_order):
                if name in bounds_rules:
                    rule = bounds_rules[name]
                    if bound_type == 'lb' and rule['lb'] is not None:
                        val = rule['lb']
                        if str(val) != bounds_values[i]:
                            bounds_updated = True
                    elif bound_type == 'ub' and rule['ub'] is not None:
                        val = rule['ub']
                        if str(val) != bounds_values[i]:
                            bounds_updated = True
                    else:
                        val = bounds_values[i]
                    new_values.append(str(val))
                else:
                    new_values.append(bounds_values[i])
            
            if bounds_updated:
                print(f"  Updated {bound_type} bounds")
            
            new_str = ', '.join(new_values)
            content = (
                content[:match.start()] +
                f"{match.group(1)}{new_str}{match.group(3)}" +
                content[match.end():]
            )
        
        return content
    
    # Update bounds if optimization parameters are detected
    param_order = detect_optimized_parameters(content)
    if param_order:
        content = update_bounds_array(content, 'lb', bounds_rules, param_order)
        content = update_bounds_array(content, 'ub', bounds_rules, param_order)
    
    # Write the updated content back to the file
    try:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"  ✓ Successfully updated {updated_count} parameters")
        return True
    except Exception as e:
        print(f"  ✗ Error writing file: {e}")
        return False

def process_folder_recursively():
    """
    Recursively processes all subfolders to find and update modified MATLAB files.
    """
    
    # Define parameter sets for each file type
    parameter_sets = {
        'blattner': {
            'r_bc': 0.019,
            'r_bp': 0.034,
            'r_cp': 0.0154,
            'sigma_bc': 1.660,
            'sigma_bp': 1.816,
            'sigma_cp': 5.740,
            'sigma_p': 3.610,
            'A': 84.523,
            'sigma_A': 0.633,
            'r_p': 0.298,
        },
        'lucey': {
            'r_bc': 0.062,
            'r_bp': 0.040,
            'r_cp': 0.0156,
            'sigma_bc': 1.002,
            'sigma_bp': 2.479,
            'sigma_cp': 4.087,
            'sigma_p': 4,
            'A': 55.780,
            'sigma_A': 0.485,
            'r_p': 0.300,
        },
        'liu': {
            'r_bc': 0.015,
            'r_bp': 0.014,
            'r_cp': 0.00320,
            'sigma_bc': 1.110,
            'sigma_bp': 1.297,
            'sigma_cp': 6.552,
            'sigma_p': 3.055,
            'A': 14.450,
            'sigma_A': 0.750,
            'r_p': 0.475,
        },
        'global': {
            'r_bc': 0.038,
            'r_bp': 0.014,
            'r_cp': 0.00537,
            'sigma_bc': 1.131,
            'sigma_bp': 1.768,
            'sigma_cp': 6.100,
            'sigma_p': 4.253,
            'A': 16.203,
            'sigma_A': 0.772,
            'r_p': 0.427,
        }
    }
    
    # Get current directory
    root_dir = os.getcwd()
    
    # File patterns to look for
    file_patterns = {
        'blattner_fit.m': 'blattner',
        'lucey_fit.m': 'lucey', 
        'liu_fit.m': 'liu',
        'global_fit.m': 'global'
    }
    
    print("MATLAB Parameter Updater - Recursive Mode")
    print("=" * 60)
    print(f"Scanning directory: {root_dir}")
    print("=" * 60)
    
    # Statistics tracking
    total_files_found = 0
    total_files_updated = 0
    folders_processed = []
    
    # Walk through all subdirectories
    for root, dirs, files in os.walk(root_dir):
        # Skip the root directory itself
        if root == root_dir:
            continue
            
        folder_name = os.path.basename(root)
        files_in_folder = []
        
        # Check for target files in current folder
        for filename, file_type in file_patterns.items():
            if filename in files:
                file_path = os.path.join(root, filename)
                files_in_folder.append((file_path, file_type, filename))
                total_files_found += 1
        
        # Process files found in this folder
        if files_in_folder:
            print(f"\n📁 Processing folder: {folder_name}")
            folders_processed.append(folder_name)
            
            for file_path, file_type, filename in files_in_folder:
                params = parameter_sets[file_type]
                success = update_matlab_parameters(file_path, params, file_type)
                if success:
                    total_files_updated += 1
    
    # Print summary
    print("\n" + "=" * 60)
    print("PROCESSING COMPLETE")
    print("=" * 60)
    print(f"Folders processed: {len(folders_processed)}")
    print(f"Files found: {total_files_found}")
    print(f"Files successfully updated: {total_files_updated}")
    
    if folders_processed:
        print(f"\nFolders processed:")
        for folder in sorted(folders_processed):
            print(f"  • {folder}")
    
    if total_files_found == 0:
        print("\nNo modified MATLAB files found in subdirectories.")
        print("Expected files: blattner_modified.m, lucey_modified.m, liu_modified.m")
    elif total_files_updated < total_files_found:
        print(f"\nWarning: {total_files_found - total_files_updated} files had errors during update.")

def main():
    """
    Main function - runs the recursive folder processing automatically.
    """
    try:
        process_folder_recursively()
    except KeyboardInterrupt:
        print("\n\nProcess interrupted by user.")
    except Exception as e:
        print(f"\nUnexpected error: {e}")

if __name__ == "__main__":
    main()
