import numpy as np
import matplotlib.pyplot as plt

error55 = np.array([0.669277,0.119006,0.011904,0.000666622,3.03275e-05,8.96797e-07])
error60 = np.array([0.612051,0.0796375,0.00824817,0.000267065,1.18368e-05,1.91843e-07])
error90 = np.array([0.345196,0.00960835,0.000365634,1.89144e-06, 1.39388e-08])

N6 = np.array([i for i in range(len(error60))])
N9 = np.array([i for i in range(len(error90))])

leng = 5
print(np.polyfit(N6, np.log(error55), 1))
print(np.polyfit(N6, np.log(error60), 1))
print(np.polyfit(N9, np.log(error90), 1))

# rate60 = [-i for i in range(len(error60))]
# rate90 = [np.cos(.9*np.pi)*rate60[i]/np.cos(.6*np.pi) for i in range(len(error60))]

plt.plot(np.log10(error60), label="0.6")
# plt.plot(rate60, label="rate 0.6")
plt.plot(np.log10(error90), label="0.9")
plt.legend()
plt.show()
