import math

data = [
    (50, 736, 157),
    (64, 721, 228),
    (55, 810, 170),
    (35, 865, 256),
    (21, 604, 443.5),
    (95, 699, 417),
    (34, 651, 178),
    (35, 1020, 230),
    (304, 675, 285.8),
]

tn = 0  # total n
tx = 0  # total sum of x
txx = 0  # total sum of x squared

for n, mean, sd in data:
    sum_x = mean * n
    sum_x2 = (sd**2 * (n-1)) + ((sum_x**2) / n)
    
    tn += n
    tx += sum_x
    txx += sum_x2

combined_mean = tx / tn
combined_sd = math.sqrt((txx - (tx**2 / tn)) / (tn - 1))

print(f"Combined Sample Size: {tn}")
print(f"Combined Mean: {combined_mean:.2f}")
print(f"Combined Standard Deviation: {combined_sd:.2f}")
print(f"Final Result: {combined_mean:.2f} ± {combined_sd:.2f} pg/ml")
