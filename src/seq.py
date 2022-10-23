import matplotlib.pyplot as plt
import numpy as np

a0 = 10

size = 5000

n = range(1,size)


a = []
b = []
c = []

for num in n:
    a_num = a0 + np.sqrt(num)
    b_num = a0 + np.log2(num**2)
    c_num = a0 + np.log10(num)
    # print(a_num)
    a.append(a_num)
    b.append(b_num)
    c.append(c_num)


# print(a)
dists_a = []
dists_b = []
dists_c = []

for i in range(len(a)-1):
    num_actual = a[i+1]
    if i >= 2:
        da = a[i+1]-a[i]
        db = b[i+1]-b[i]
        dc = c[i+1]-c[i]
        # print(f"La distancia entre el término {i+1} y el término {i} es ", d)
        dists_a.append(da)
        dists_b.append(db)
        dists_c.append(dc)

# print(len(dists))

x = range(1, size-3)

fig, ax = plt.subplots()

ax.plot(x, dists_a, label='Distancia entre términos de A')
ax.plot(x, dists_b, label='Distancia entre términos de B')
ax.plot(x, dists_c, label='Distancia entre términos de C')
ax.legend()

plt.show()
