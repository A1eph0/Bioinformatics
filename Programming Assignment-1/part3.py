import matplotlib.pyplot as plt
import parameters # importing necessary macro values

genome_unedited = input("Enter genome sequence:\n")
genome = genome_unedited.replace("\n", "")

# making the restriction map
restriction_map = [i+1 for i in range(len(genome)) if genome.startswith(parameters.ecor1_sequence, i)]

print("Restriction map for EcoR1 is as given below:", restriction_map, sep="\n", end="\n")


# plotting the restriction sites on restriction map
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(0,len(genome))
ax.set_ylim(0,7)

xmin = 0
xmax = len(genome)
y = 5
height = 1

plt.title("Restriction map with EcoR1")

plt.hlines(y, xmin, xmax)

plt.vlines(xmin, y - height / 2., y + height / 2.)
plt.vlines(xmax, y - height / 2., y + height / 2.)
plt.text(xmin, y, '0', horizontalalignment='right')
plt.text(xmax, y, str(len(genome)), horizontalalignment='left')

flag = 1
for i in restriction_map:
    plt.vlines(i, y - height / 2., y + height / 2.)
    plt.text(i, y - flag*((height / 2.) + 0.5) , str(i), horizontalalignment='center')
    flag*=-1

plt.axis('off')
plt.show()