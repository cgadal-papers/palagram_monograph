import glob
import os

list_script = sorted(glob.glob('figure*.py', recursive=True))

for script in list_script:
    print(script)
    os.system('python3 {}'.format(script))
