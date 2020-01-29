import os, platform, sys

OS = platform.system()
Version = platform.version().split()[0]
PyVersion = sys.version.split()[0]
home = os.getcwd()

print('Runing Setup for MoleKing_util on {} -{}- with python {}:'.format(OS, Version, PyVersion))

try:
    os.chdir('{}/MoleKing_util'.format(home))
    for arq in os.listdir():
        os.remove(arq)
    os.chdir(home)
    os.removedirs('MoleKing_util')
except:
    pass

d = '{}/src'.format(home)
CList = ['{}/main.cpp'.format(d)]
subdirs = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]
for directory in subdirs:
    for arq in os.listdir(directory):
        if arq.split('.')[1] == 'cpp':
            CList.append('{}/{}'.format(directory, arq))

os.makedirs('MoleKing_util')

if OS == 'Linux':
    flags = '-O3 -Wall -shared -std=c++11 -fPIC `python{} -m pybind11 --includes`'.format(PyVersion[0:3])
    target = 'MoleKing_util`python{}-config --extension-suffix`'.format(PyVersion[0:3])
elif OS == 'Darwin':
    flags = '-O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes`'
    target = 'MoleKing_util`python3-config --extension-suffix`'
else:
    print('Recomendamos usar um OS feito para adultos...')

objs = ' '.join(CList)
os.chdir('{}/MoleKing_util'.format(home))
obj = os.popen('c++ {0} {2} -o {1}'.format(flags, target, objs), 'r')
obj.read()
print('Success, MoleKiing_util was compiled on {} -{}- with python {}!'.format(OS, Version, PyVersion))







###### Comments ######
'''for Cfile in CList:
    if 'main' in Cfile:
        obj = os.popen('c++ -o {0}.o {1}.cpp {2}'.format(Cfile.split('/')[-1], Cfile, flags), 'r')
        obj.read()
    else:
        obj = os.popen('c++ -o {0}.o {1}.cpp {1}.hh {2}'.format(Cfile.split('/')[-1], Cfile, flags), 'r')
        obj.read()'''