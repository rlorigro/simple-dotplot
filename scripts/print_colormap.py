from matplotlib import colors
from matplotlib import cm
from matplotlib import pyplot
import argparse
import numpy


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap


def main():
    colormap = cm.get_cmap('coolwarm', 256)

    a = pyplot.axes()
    colors = list()
    x = list()
    y = list()

    for i in range(256):
        c = float(i)/256
        print("{%f,%f,%f}," % (colormap(c)[0],colormap(c)[1],colormap(c)[2]))

        colors.append(colormap(c))
        x.append(i)
        y.append(1)

    pyplot.scatter(x,y,color=colors)
    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()
