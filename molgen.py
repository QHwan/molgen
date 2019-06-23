import copy
import math
import numpy as np

class Atom:
    def __init__(self, name, position):
        self.atom_name = name
        self.atom_position = position


class Residue:
    def __init__(self, name):
        self.residue_name = name
        self.atom_vec = []

    def add_atom(self, atom):
        self.atom_vec.append(atom)

    def make_periodic(self, nx=1, ny=1, nz=1, dx=0, dy=0, dz=0):
        ref_atom_vec = copy.deepcopy(self.atom_vec)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for ref_atom in ref_atom_vec:
                        atom_name = ref_atom.atom_name
                        ref_x, ref_y, ref_z = ref_atom.atom_position
                        x = ref_x + i*dx
                        y = ref_y + j*dy
                        z = ref_z + k*dz
                        atom = Atom(atom_name, [x, y, z])
                        self.atom_vec.append(atom)


class Total:
    def __init__(self, box, title="default system"):
        self.title = title
        self.box = box
        self.residue_vec = []

    def add_residue(self, residue):
        self.residue_vec.append(residue)

    def get_num_atoms(self):
        num_atoms = 0
        for residue in self.residue_vec:
            num_atoms += len(residue.atom_vec)
        return num_atoms


def write_gro(ofilename, total):
    with open(ofilename, "w") as ofile:
        ofile.write(str(total.title) + '\n')
        ofile.write(str(total.get_num_atoms()) + '\n')
        res_num = 0
        res_name = ""
        atom_name = ""
        atom_num = 0
        x = 0
        y = 0
        z = 0
        for residue in total.residue_vec:
            res_num += 1
            res_name = residue.residue_name
            for atom in residue.atom_vec:
                atom_name = atom.atom_name
                atom_num += 1
                x, y, z = atom.atom_position
                ofile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (res_num, res_name, atom_name, atom_num, x, y, z))
                ofile.write("\n")

        ofile.write("%.3f\t%.3f\t%.3f" % (box[0], box[1], box[2]))
        ofile.write("\n")

def make_graphene(total, x, y, set_box=True):
    GRA = Residue("GRA")

    ccbond = 0.142
    dx=ccbond*math.cos(30*(math.pi/180.0))
    dy=ccbond*math.sin(30*(math.pi/180.0))  
    C1 = Atom("C", [0., 0., 0.])
    C2 = Atom("C", [-1*dx, dy, 0.])
    C3 = Atom("C", [-1*dx, dy+ccbond, 0.])
    C4 = Atom("C", [0, 2*dy+ccbond, 0.])
    GRA.add_atom(C1)
    GRA.add_atom(C2)
    GRA.add_atom(C3)
    GRA.add_atom(C4)
    nx = int(x/(2*dx))
    ny = int(y/(2*dy+2*ccbond))
    GRA.make_periodic(nx=nx, ny=ny, nz=1, dx=2*dx, dy=2*dy+2*ccbond, dz=0)
    tot.add_residue(GRA)

    if set_box == True:
        tot.box[0] = nx*2*dx
        tot.box[1] = ny*2*(dy+ccbond)


box = [10, 10, 10]
tot = Total(box)
make_graphene(tot, 5, 5)
ofilename = "test.gro"
write_gro(ofilename, tot)
