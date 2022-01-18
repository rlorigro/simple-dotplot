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


def main(input_path):
    matrix = numpy.loadtxt(input_path, delimiter=',', dtype=float)
    matrix = numpy.log10(matrix+1)

    pyplot.imshow(matrix, cmap=truncate_colormap(cm.get_cmap('viridis', 256), minval=0.35, maxval=1.0))
    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        type=str,
        required=True,
        help="Path of PAF file"
    )

    args = parser.parse_args()

    main(args.i)
