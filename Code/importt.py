import os
#os.system("sudo apt update")
#os.system("sudo apt install python3-pip")
os.system("pip install openpyxl")
try:
    import Bio
except ImportError:
    os.system("pip install biopython")
    import Bio


try:
    import mne
except ImportError:
    os.system("pip install mne")
    import mne



os.system("pip3 install pandas")
try:
    import panel as pn
except:
    os.system("pip install panel")
    import panel as pn
try:
    import openpyxl
except ImportError:
    os.system("pip3 install openpyxl")


try:
    import pandas as pd
except ImportError:
    os.system("pip3 install pandas")
    import pandas as pd

try:
    import glob
except:
    os.system("sudo pip3 install glob3")
    import glob
try:
    import shutil
except:
    os.system("pip3 install shutil")
    import shutil
try:
    import glob
except:
    os.system("pip3 install glob")
    import glob
try:
    import urllib.request
except:
    os.system("pip3 install urllib.request")
try:
    import sys
except ImportError:
    os.system("pip3 install sys")
try:
    import time
except ImportError:
    os.system("pip3 install time")

try:
    import Bio
except ImportError:
    os.system("pip3 install Bio")

try:
    import subprocess
except ImportError:
    os.system("pip3 install subprocess")

try:
    import threading
except ImportError:
    os.system("pip3 install threading")

try:
    import pexpect
except ImportError:
    os.system("pip3 install pexpect")

try:
    import os.path
except ImportError:
    os.system("pip3 install os.path")

try:
    import shutil
except ImportError:
    os.system("pip3 install shutil")
try:
    import scipy.stats
except ImportError:
    os.system("pip3 install scipy matplotlib")
try:
    import scipy.stats 
except ImportError:
    os.system("python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose")
 
try:
    import matplotlib.pyplot as plt
except ImportError:
    os.system('pip install matplotlib')
    import matplotlib.pyplot as plt
import os
import time

import subprocess
import threading
import pexpect
import os.path
import shutil
import scipy.stats
from Bio import Entrez
from Bio.Seq import Seq
from subprocess import Popen, PIPE
from Bio.Phylo.PAML import codeml
from Bio import AlignIO
import pandas as pd
from scipy.stats import chi2
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os, io, random
import string
import numpy as np

from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
   
import panel.widgets as pnw
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

