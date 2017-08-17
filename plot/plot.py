import argparse
import json
import matplotlib.pyplot as plt

def _main():
    parser = argparse.ArgumentParser("Plot total DOS")
    parser.add_argument("prefix", type=str, help="Calculation prefix")

    args = parser.parse_args()

    with open("{}_dos.json".format(args.prefix)) as fp:
        dos_data = json.load(fp)

    es = dos_data["es"]
    total_dos = dos_data["total_dos"]

    plt.plot(es, total_dos, 'k-')
    plt.show()

if __name__ == "__main__":
    _main()
