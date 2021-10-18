import matplotlib.pyplot as plt

# setting alpha value
alpha = 0.005

# initialisting start conditions
seq_len = 1000
num_correct = seq_len
num_incorrect = 0

# array for storing plot points
data = [0]

# looping over each time step
for _ in range(500):
    new_incorrect = int( 3 * alpha * num_correct)
    new_correct = int( alpha * num_incorrect)

    num_correct = num_correct - new_incorrect + new_correct
    num_incorrect = seq_len - num_correct

    data.append(num_incorrect/seq_len)


# plot graph
plt.plot(data)
plt.show()
