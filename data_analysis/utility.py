from matplotlib import pyplot as plt

def plot():
    fig, ax = plt.subplots()
    ax.grid()
    ax.tick_params(which='both', direction="in")

    return fig, ax

def plot_scattering(ax, x, y_re, y_im, label=""):
    lines = ax.plot(x, y_re, label=label)
    ax.plot(x, y_im, linestyle="--", color=lines[0].get_color())