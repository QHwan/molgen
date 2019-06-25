import copy
import math
from random import uniform
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
                    if i == 0 and j == 0 and k == 0:
                        continue
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

    def make_replicate(self, residue, nx=1, ny=1, nz=1, dx=0, dy=0, dz=0):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    new_residue = copy.deepcopy(residue)
                    for atom in new_residue.atom_vec:
                        atom.atom_position[0] = ref_x + i*dx
                        atom.atom_position[1] = ref_y + j*dy
                        atom.atom_position[2] = ref_z + k*dz
                    self.residue_vec.append(new_residue)

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

        ofile.write("%.3f\t%.3f\t%.3f" % (total.box[0], total.box[1], total.box[2]))
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
    total.add_residue(GRA)

    if set_box == True:
        total.box[0] = nx*2*dx
        total.box[1] = ny*2*(dy+ccbond)

def make_cristobalite(total, x, y, set_box=True, bind_hydrogen=False):
    CRI = Residue("CRI")
  
    SI1 = Atom("SI", [0., 0., 0.])
    #SH1 = Atom("SH", [-0.246, 0.142, -0.1])
    SI5 = Atom("SI", [-0.246, 0.142, -0.1])
    SI2 = Atom("SI", [-0.246, 0.426, 0.])
    #SH2 = Atom("SH", [0., 0.568, -0.1])
    SI6 = Atom("SI", [0., 0.568, -0.1])
    O1 = Atom("O", [-0.123, 0.071, -0.05])
    O2 = Atom("O", [0., 0., 0.151])
    O3 = Atom("O", [0.123, 0.071, -0.05])
    O4 = Atom("O", [-0.246, 0.284, -0.05])
    O5 = Atom("O", [0.123, 0.497, -0.05])
    O6 = Atom("O", [-0.246, 0.426, 0.151])
    O7 = Atom("O", [-0.123, 0.497, -0.05])
    O8 = Atom("O", [0, 0.71, -0.05])
    SI3 = Atom("SI", [0., 0., 0.302])
    SI4 = Atom("SI", [-0.246, 0.426, 0.302])
    SH3 = Atom("SH", [-0.246, 0.142, 0.402])
    SH4 = Atom("SH", [0., 0.568, 0.402])
    O9 = Atom("O", [-0.123, 0.071, 0.352])
    O10 = Atom("O", [0.123, 0.071, 0.352])
    O11 = Atom("O", [-0.246, 0.284, 0.352])
    O12 = Atom("O", [0.123, 0.497, 0.352])
    O13 = Atom("O", [-0.123, 0.497, 0.352])
    O14 = Atom("O", [0., 0.71, 0.352])
    #OH1 = Atom("OH", [-0.246, 0.142, -0.251])
    #OH2 = Atom("OH", [-0, 0.567, -0.251])
    OH1 = Atom("OH", [-0.246, 0.142, 0.553])
    OH2 = Atom("OH", [0., 0.568, 0.553])
    CRI.add_atom(SI1)
    #CRI.add_atom(SH1)
    CRI.add_atom(SI2)
    #CRI.add_atom(SH2)
    CRI.add_atom(SI5)
    CRI.add_atom(SI6)
    CRI.add_atom(O1)
    CRI.add_atom(O2)
    CRI.add_atom(O3)
    CRI.add_atom(O4)
    CRI.add_atom(O5)
    CRI.add_atom(O6)
    CRI.add_atom(O7)
    CRI.add_atom(O8)
    CRI.add_atom(SI3)
    CRI.add_atom(SI4)
    CRI.add_atom(SH3)
    CRI.add_atom(SH4)
    CRI.add_atom(O9)
    CRI.add_atom(O10)
    CRI.add_atom(O11)
    CRI.add_atom(O12)
    CRI.add_atom(O13)
    CRI.add_atom(O14)
    #CRI.add_atom(OH1)
    #CRI.add_atom(OH2)
    CRI.add_atom(OH1)
    CRI.add_atom(OH2)

    if bind_hydrogen == True:
        theta = uniform(0,360)
        dx = 0.1*math.cos(theta*3.14/180)
        dy = 0.1*math.sin(theta*3.14/180)
        HH1 = Atom("HH", [-0.246+dx, 0.142+dy, 0.586])
        HH2 = Atom("HH", [0.+dx, 0.568+dy, 0.586])
        CRI.add_atom(HH1)
        CRI.add_atom(HH2)

    dx = 0.492
    dy = 0.852
    nx = int(x/dx)
    ny = int(y/dy)
    print(nx, ny)
    CRI.make_periodic(nx=nx, ny=ny, nz=1, dx=dx, dy=dy, dz=0)
    total.add_residue(CRI)

    if set_box == True:
        total.box[0] = nx*dx
        total.box[1] = ny*dy

def make_water(total, position, model="SPC/E"):
    SOL = Residue("SOL")
    x, y, z = position

    if model == "SPC/E":
        OW = Atom("OW", [0. + x, 0. + y, 0. + z])
        HW1 = Atom("HW1", [0.099 + x, -0.011 + y, -0.009 + z])
        HW2 = Atom("HW2", [-0.04 + x, -0.086 + y, 0.032 + z])
    SOL.add_atom(OW)
    SOL.add_atom(HW1)
    SOL.add_atom(HW2)

    total.add_residue(SOL)



box = [10, 10, 10]
tot = Total(box)
make_cristobalite(tot, 5, 5.2, bind_hydrogen=True)

water_nx = int(4.5/0.3)
water_ny = int(4.5/0.3)
water_nz = int(5./0.3)
for i in range(water_nx):
    for j in range(water_ny):
        for k in range(water_nz):
            make_water(tot, [0.3*i, 0.3*j, 0.8+0.3*k], model="SPC/E")

ofilename = "cri_sol_q1.0.gro"
write_gro(ofilename, tot)
