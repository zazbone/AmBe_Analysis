# Usage: $ python generate_utils.py input/file/path.root [output/file/path.root (Optional)]


import ROOT
from _path import *
from utils import load_AmBe
from sys import argv
import os

load_AmBe()

if len(argv) < 2:
    # TODO: Find the write exception to raise
    raise ("input files not provided, abort")

input_file = Path(argv[1])
if len(argv) == 3:
    output_file = Path(argv[2])
else:
    print("Output path not provided, tmp folder created")
    os.mkdir("tmp/")
    output_file = Path("./tmp") / f"{input_file.name}"

ROOT.EQuenched(str(input_file), str(output_file))