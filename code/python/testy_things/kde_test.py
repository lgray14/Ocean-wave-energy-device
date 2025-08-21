import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

# Generate random data
data = np.random.normal(0, 1, size=20)

# Create a Gaussian KDE object
kde = gaussian_kde(data)
# print(kde)

# Evaluate the density over a range
x = np.linspace(-5, 5, 1000)
density = kde(x)
# print(density)

# Plot the results
plt.plot(x, density, label="KDE")
plt.hist(data, bins=30, density=True, alpha=0.5, label="Histogram")
plt.legend()
plt.title("Kernel Density Estimation")
plt.show()