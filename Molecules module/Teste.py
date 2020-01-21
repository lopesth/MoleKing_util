import molecules
x = molecules.Molecule()

x.addAtom('C',-0.73876630,2.99880706,1.65325738)
x.addAtom('H',-0.80221082,2.45489630,0.70650216)
x.addAtom('H',-1.42674195,2.51595206,2.35307903)
x.addAtom('H',-1.10628144,4.01243982,1.47627677)
x.addAtom('C', 0.68662235,3.00063539,2.19383677)
x.addAtom('H', 1.34401523,3.51973092,1.48863864)
x.addAtom('H', 0.72448594,3.58022076,3.12205084)
x.addAtom('C', 1.22222258,1.59635554,2.44896950)
x.addAtom('H', 2.24394207,1.61645041,2.83565590)
x.addAtom('H', 0.60429125,1.06343748,3.17723225)
x.addAtom('H', 1.22868677,1.00240854,1.53060295)

'''y = molecules.Atom('C',-0.73876630,2.99880706,1.65325738)
print(y.getX())
y.setX(1)
print(y.getX())



for atom in x.getMolecule():
    print('{}        {}    {}    {}'.format(*atom))

x.stdOrientation()
print('######################################\n')
'''
x.moveMassCenter()
for atom in x.getMolecule():
    print('{}        {}    {}    {}'.format(*atom))
