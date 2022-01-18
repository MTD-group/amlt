import numpy as np
from ase import Atoms, io


from amlt.kgrid import get_kpts_from_kpd

def demo_structure(
            La,
            Lb,
            Lc,
            angle1,
            angle2):

    cell = [[                         La,                           0,                           0],
            [Lb*np.cos(angle1*np.pi/180), Lb*np.sin(angle1*np.pi/180),                           0],
            [Lc*np.cos(angle2*np.pi/180),                           0, Lc*np.sin(angle2*np.pi/180)]]
    atoms = Atoms('H', cell=cell)

    return atoms
    
def print_info(atoms):
    rcell = atoms.cell.reciprocal()
    l  = np.linalg.norm(atoms.cell,axis=1)
    rl = np.linalg.norm(atoms.cell.reciprocal(),axis=1)
    print('cell vol', atoms.get_volume())
    print('cell\n',   l[0], atoms.cell[0], '\n',  l[1], atoms.cell[1], '\n',  l[2], atoms.cell[2])
    print('rcell\n', rl[0],      rcell[0], '\n', rl[1],      rcell[1], '\n', rl[2],      rcell[2])


############
print('-----test 1-----')
kpd = 8000
atoms = demo_structure(            
            La    =  2,
            Lb    =  4,
            Lc    = 100,
            angle1=30,
            angle2=90)
io.write('kgrid_test_atoms_1.CONTCAR', atoms)
print_info(atoms)
kpts = get_kpts_from_kpd(atoms,kpd)
print('-----test 1 even-----')
kpts = get_kpts_from_kpd(atoms,kpd,only_even=True)
print()



#############
print('-----test 2-----')
kpd = 9000
atoms = demo_structure( 
            La =  3,
            Lb =  6,
            Lc = 1,
            angle1=30,
            angle2=90)
            
io.write('kgrid_test_atoms_2.CONTCAR', atoms)
print_info(atoms)
kpts = get_kpts_from_kpd(atoms,kpd)
print('-----test 2 even-----')
kpts = get_kpts_from_kpd(atoms,kpd,only_even=True)
print()





#############
print('-----test 3-----')
kpd = 4096
atoms = demo_structure( 
            La =  3,
            Lb =  3,
            Lc = 64,
            angle1=90,
            angle2=90)
            
io.write('kgrid_test_atoms_3.CONTCAR', atoms)
print_info(atoms)
kpts = get_kpts_from_kpd(atoms,kpd)
print('-----test 3 even-----')
kpts = get_kpts_from_kpd(atoms,kpd,only_even=True)
print()
