try:
    print(_PATH_SET)
except NameError:
    import os as _os
    from pathlib import Path as _Path
    _os.chdir(_Path(__file__).parent / "../../script")
    import ROOT
    #import utils
