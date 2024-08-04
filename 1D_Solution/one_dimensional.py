import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt

def init_conditions(h):
    """ Set Initial Conditions to a parabola """
    init = np.arange(-1, 1, h, dtype=float)
    x = np.zeros((100, init.shape[0]))
    x[0, :] = init
    x[0, :] = 0.5 - x[0, :]**2
    return x

def direction_matrix(x):
    """ Unused """
    return np.tile(np.array(-3), x.shape[0])

def plotter(x, h):
    """ Plot the function """
    def animate(i):
        line.set_ydata(data[i, :])  # update the data.
        return line

    data = x.copy()
    fig, ax = plt.subplots()
    x = np.arange(-1, 1, h, dtype=float)
    line, = ax.plot(x, data[0, :])
    ax.set_ylim([0, 1])
    ax.set_title('1D Lax-Friedrichs Method to Non-Linear PDE')
    ax.set_xlabel('Space Discretization')
    ax.set_ylabel('Time Discretization')


    ani = animation.FuncAnimation(
        fig, animate, interval=300)

    ani.save('1D_plot.mp4', fps=10)
    #plt.show()


def main():
    k = 0.01
    h = 0.1

    x = init_conditions(h)

    # CFL Condition Check
    if abs((2 * -3 * k) / h) > 1:
        raise ValueError("Sim Unstable")

    # Iterate over each time bin calculating Lax-Friedrichs Solution
    for t in range(0, 99):
        av = 0.5 * (x[t, 2:] + x[t, :-2])
        lax = (x[t, 2:] * 3 * (1 - x[t, 2:]/1) - x[t, :-2] * 3 * (1 - x[t, :-2]/1))
        x[t+1, 1:-1] = av - (0.01 / (2 * h)) * lax
        x[t+1, 0] = x[t+1, 1]
        x[t+1, -1] = x[t+1, -2]

        x[t+1, x[t+1, :] > 1] = 1

    plotter(x, h)

if __name__ == "__main__":
    main()