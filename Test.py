
import pynapple as nap
from pathlib import Path

# load NWB file

datapath = Path('/Users/doctordu/Dropbox (Personal)/Professional/Data/Trisha')
foldername = 'H7101-241201'

filepath = datapath / foldername
filename = foldername+ '.nwb'
filepath = filepath / filename

data = nap.load_file(str(filepath))
spikes = data['units']
epochs = data['epochs']
ang = data['head-direction']
pos = data['position']


print()