import importlib.machinery
import types
import pandas as pd
import os

loader = importlib.machinery.SourceFileLoader('color','/Users/bgprocess/aroma/AromaPack/finder_colors.py')
color = types.ModuleType(loader.name)
loader.exec_module(color)
# color.get('/Users/bgprocess/aroma/AromaPack/finder_colors.py')
# color.set('/Users/bgprocess/aroma/AromaPack/finder_colors.py','red')

workingdir = '/Users/bgprocess/aroma/hg19'
dir = os.path.join(workingdir,'processed')
files = os.listdir(dir)
files = [f for f in files if f.startswith('aroma')]
print(files)
ls = []
for f in files:
    if os.stat(os.path.join(dir,f)).st_size != 0:
        print (f)
        aroma = pd.read_table(os.path.join(dir,f),header=None)
        aroma.shape

        e = aroma.iloc[:,4]
        for i,p in enumerate(e):
            if p == 'Error in if (sampleAA * sampleBB > 0L) {: argument is of length zero' or p.startswith('Error: [affxparser Fusion SDK exception]'):
                ls.append(aroma.iloc[i,3])
ls = list(set(ls))
print(ls)
for s in ls:
    color.set('/Volumes/arraymapMirror/arraymap/hg19/{0}'.format(s),'gray')
