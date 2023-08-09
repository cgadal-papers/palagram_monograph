import glob
import os

list_script = sorted(glob.glob('figure*.py', recursive=True))

for script in list_script[1:]:
    print(script)
    os.system('python3 {}'.format(script))
