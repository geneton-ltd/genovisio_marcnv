import argparse
import json
import os
import sys
from dataclasses import asdict

import annotation

from marcnv.src.acmg.acmg_classify import MarCNVClassifier


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify CNV and/or find intersecting items in MongoDB collections.")
    parser.add_argument("input", help="Annotated CNV stored as json")
    parser.add_argument("--output", help="Path to store the prediction JSON. Else prints to stdout.", default=None)
    args = parser.parse_args()

    whole_cnv_annotation = annotation.Annotation.load_from_json(args.input)
    classifier = MarCNVClassifier(whole_cnv_annotation)
    prediction = classifier.classify()

    if args.output:
        path = os.path.abspath(args.output)
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump(asdict(prediction), f, indent=2)
    else:
        print(json.dumps(asdict(prediction), indent=2), file=sys.stdout)


if __name__ == "__main__":
    main()
