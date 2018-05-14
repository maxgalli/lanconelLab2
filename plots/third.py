import matplotlib.pyplot as plt

contrast = [0.3, 0.38, 0.5, 0.6]
error = [0.1, 0.09, 0.08, 0.07]
thickness = [4, 6, 8, 10]

plt.xlabel("Box thickness (mm)")
plt.ylabel("Contrast")

plt.errorbar(thickness, contrast, yerr=error, marker='o')

plt.savefig("contrast.png")