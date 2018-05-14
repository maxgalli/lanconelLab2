import matplotlib.pyplot as plt

# Look at the ods for reference (Because I am an Horrible Person)
data = [[0.99831,0.70302,0.63643,0.08849],
[0.99853,0.91125,0.81241,0.18645],
[0.99822,0.98782,0.85712,0.35822],
[0.99848,0.99394,0.85964,0.68205],
[0.99843,0.99414,0.85900,0.89394]]

tdata = []
for j in range(0,4):
	temp = []
	for i in range(0, 5):
		temp.append(data[i][j])
	tdata.append(temp)

print(tdata)

kevs = [10, 30, 50, 140]
thickness = [0.5, 1, 2, 5, 10]

#%%

plt.ylabel("QDE")
plt.xlabel("Photon energy (KeV)")

plt.plot(kevs, data[0][::1], marker = "o", label = "0.5 mm")
plt.plot(kevs, data[1][::1], marker = "o", label = "1 mm")
plt.plot(kevs, data[2][::1], marker = "o", label = "2 mm")
plt.plot(kevs, data[3][::1], marker = "o", label = "5 mm")
plt.plot(kevs, data[4][::1], marker = "o", label = "10 mm")

plt.legend(title = "Detector thickness:",loc=3)
plt.savefig("qde_en.png")

plt.clf()

#%%

plt.ylabel("QDE")
plt.xlabel("Detector thickness (mm)")

plt.plot(thickness, tdata[0], marker = "o", label = "10 Kev")
plt.plot(thickness, tdata[1], marker = "o", label = "30 Kev")
plt.plot(thickness, tdata[2], marker = "o", label = "50 Kev")
plt.plot(thickness, tdata[3], marker = "o", label = "140 Kev")

plt.ylim(0,1.1)

plt.legend(title = "Photon energy:",loc=4)
plt.savefig("qde_th.png")

plt.clf()
