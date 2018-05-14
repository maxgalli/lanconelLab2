import matplotlib.pyplot as plt

# Look at the ods for reference (Because I am an Horrible Person)
data = [[0.000056095, 0.000088196, 0.006918400, 0.019660100],
[0.037553000, 0.024839500, 0.032030500, 0.029375000],
[0.000068127, 0.003997350, 0.007198050, 0.021206400],
[0.076704000, 0.056588000, 0.064437000, 0.064977500],
[0.000133245, 0.004027900, 0.007315550, 0.022757400],
[0.154841500, 0.129461500, 0.136605500, 0.090522000],
[0.000133598, 0.001587190, 0.007343750, 0.058068500],
[0.389630000, 0.361665000, 0.367775000, 0.250510000],
[0.000133598, 0.008213250, 0.007423650, 0.060606500],
[0.780905000, 0.752235000, 0.757640000, 0.560240000]]

tdata = []
for j in range(0,4):
	temp = []
	for i in range(0, 10):
		temp.append(data[i][j])
	tdata.append(temp)

print(tdata)

kevs = [10, 30, 50, 140]
thickness = [0.5, 1, 2, 5, 10]

#%%

plt.ylabel("FWHM (cm)")
plt.xlabel("Photon energy (KeV)")

plt.plot(kevs, data[1], marker = "o", label = "0.5 mm")
plt.plot(kevs, data[3], marker = "o", label = "1 mm")
plt.plot(kevs, data[5], marker = "o", label = "2 mm")
plt.plot(kevs, data[7], marker = "o", label = "5 mm")
plt.plot(kevs, data[9], marker = "o", label = "10 mm")

plt.legend(title = "Detector thickness:")
plt.savefig("kev_fwhm_blur.png")

plt.clf()

#%%

plt.ylabel("FWHM (cm)")
plt.xlabel("Photon energy (KeV)")

plt.plot(kevs, data[0], marker = "o", label = "0.5 mm")
plt.plot(kevs, data[2], marker = "o", label = "1 mm")
plt.plot(kevs, data[4], marker = "o", label = "2 mm")
plt.plot(kevs, data[6], marker = "o", label = "5 mm")
plt.plot(kevs, data[8], marker = "o", label = "10 mm")

plt.legend(title = "Detector thickness:")
plt.savefig("kev_fwhm_noblur.png")

plt.clf()

#%%

plt.ylabel("FWHM (cm)")
plt.xlabel("Detector thickness (mm)")

plt.plot(thickness, tdata[0][::2], marker = "o", label = "10 Kev")
plt.plot(thickness, tdata[1][::2], marker = "o", label = "30 Kev")
plt.plot(thickness, tdata[2][::2], marker = "o", label = "50 Kev")
plt.plot(thickness, tdata[3][::2], marker = "o", label = "140 Kev")

plt.legend(title = "Photon energy:")
plt.savefig("mm_fwhm_noblur.png")

plt.clf()

#%%

plt.ylabel("FWHM (cm)")
plt.xlabel("Detector thickness (mm)")

plt.plot(thickness, tdata[0][1::2], marker = "o", label = "10 Kev")
plt.plot(thickness, tdata[1][1::2], marker = "o", label = "30 Kev")
plt.plot(thickness, tdata[2][1::2], marker = "o", label = "50 Kev")
plt.plot(thickness, tdata[3][1::2], marker = "o", label = "140 Kev")

plt.legend(title = "Photon energy:")
plt.savefig("mm_fwhm_blur.png")

plt.clf()
