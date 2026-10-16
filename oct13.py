from random import random
import matplotlib.pyplot as plt
import numpy as np

total = 1000
inside = 0
xs = []
ys = []

for i in range(total):
    x = random()
    y = random()
    xs.append(x)
    ys.append(y)
    if (x**2+y**2) < 1:
        inside += 1
        
ratio = inside/(total * np.pi)

print(ratio)

x = np.linspace(0,1,100)
y = np.sqrt(1-x**2)

plt.figure(figsize=(6,6))
plt.plot(x,y)
plt.scatter(xs,ys,s=5,marker='^',facecolors = 'none', edgecolors='red',linewidths=0.5)

plt.xlim(0,1)
plt.ylim(0,1)



plt.show()


