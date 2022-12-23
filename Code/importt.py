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
try:
    import panel as pn
except:
    os.system("pip install panel")
    import panel as pn


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
    import os.path
except ImportError:
    os.system("pip3 install os.path")

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

