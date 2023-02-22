import argparse
from collections import defaultdict

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def parse_args():
    parser = argparse.ArgumentParser(description="Visualize Jaccard by sample")
    parser.add_argument(
        "--input",
        required=True,
        nargs="+",
        dest="files",
        help="Input files",
    )
    parser.add_argument(
        "--output",
        required=True,
        dest="output",
        help="Output file",
    )
    return parser.parse_args()


def sort_by_sample(x: tuple[str, dict[int, float]]) -> float:
    """
    Sort comparator that sorts by max jaccard similarity in a segment
    """
    return max(map(float, list(x[1].values())))

def plot_jaccard_by_sample(files, output):

    # dict -> dict -> list
    # jaccard similarity by sample and segment
    jacc = defaultdict(dict)

    with PdfPages(output) as pdf:
        for file in files:
            segment = int(file.split(".")[-2])
            with open(file, "r") as f:
                for line in f:
                    sample, jaccard = line.strip().split("\t")
                    jacc[sample][segment] = float(jaccard)

        for sample, seg_dict in sorted(jacc.items(), key=sort_by_sample, reverse=True):
            plt.figure(figsize=(10, 5))
            sns.barplot(
                x=list(seg_dict.keys()),
                y=list(seg_dict.values()),
            )
            plt.title(sample)
            plt.xlabel("Segment")
            plt.ylabel("Jaccard similarity")
            pdf.savefig()
            plt.close()

if __name__ == "__main__":
    args = parse_args()
    plot_jaccard_by_sample(args.files, args.output)
