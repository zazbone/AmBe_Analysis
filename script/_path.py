from pathlib import Path

FILE_PATH = Path(__file__).absolute()
PROJECT_PATH = FILE_PATH.parents[1]
INCLUDE_PATH = PROJECT_PATH / "include"
DLL_PATH = PROJECT_PATH / "libAmBe.so"
DATASET_PATH = PROJECT_PATH / "dataset"
PLOTS_PATH = PROJECT_PATH / "plots"
SCRIPT_PATH = PROJECT_PATH / "script"