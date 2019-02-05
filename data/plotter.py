import matplotlib.pyplot as plt


fig, ax = plt.subplots()

xs = [12, 10, 8, 6]
val = [14.385806926284722, 1.5719116803128077, 72.24672982833195, 439.40974789609237]

plt.bar(xs, val)
plt.show()