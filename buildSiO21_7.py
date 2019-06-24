#!/usr/bin/python
'''
    
    buildSiO21_3
    
    Author: KimQHwan
    Date: November 8th 2012
    ********************************************
    Description:
    Build SiO2 crystal. Output file is saved in .gro format. 
    
    
    Syntax:
    
    buildWatercluster.py [options] outfile
    
    The available options are:
    --version
    
    shows prOram's version number and exit.
    
    -c, --credits
    
    display credits.
    
    -b, --box

	specify the size of box vector(A).

    -m, --model
    
    choose the water model. Valid structures are: a-Quartz, cristobalite    
    
	-p, --pbc

	determine periodicity of flat gold. only xy periodicity is supported.
    
    Requirements:
    Require numpy installed
  
    
    + v 1.0 - July 2012:
    - first release of buildSiO2
	+ v 1.1 - August 2012:
	- periodic option is added.
	+ v 1.2 - August 2012:
	- transformation from the hexagonal unit cell -> orthorhombic unit cell
	+ v 1.3 - September 2012:
	- cristobalite structure is added - only for xy period.
    
    License:
    Freeware. You can use, modify and redistribute the source.
    
    
    Contacts:
    You liked this prOram? You have suggestions?
    Drop me a mail at kim.qhwan@gmail.com
    '''
    
#import modules
from sys import argv, exit, path, version
from optparse import OptionParser as OP 
from numpy import zeros, pi, sin, cos, tan, modf, ceil, sqrt, matrix, array
from os import path as path_os
from os import system
from random import uniform

'''
#===============================================================================
#                               SUBROUTINES
#===============================================================================
'''


def getdist(at1,at2):
	''' Calculate distance between two particles
	'''
	dist_at=sqrt((at2[0]-at1[0])**2+(at2[1]-at1[1])**2+(at2[2]-at1[2])**2)
	return dist_at

def filecheck(file):
	''' Check if infile exists
	'''
	if path_os.isfile(file) == 0:
		found=False
	else:
		found=True
	return found

def parsecmd():
	description="Build the SiO2 crystal."
	usage = "usage: %prO [options] output_file"
	#parse command line
	parser=OP(version='%prO 1.0',description=description, usage=usage)
    
	parser.add_option('-c','--credits',dest='credits',action='store_true',
					default=False,help='display credits')
	parser.add_option('-b','--box',dest='box',nargs=3,type='float',\
					help='determine the length of box vector')
	parser.add_option('-m','--model',dest='model',nargs=1,type='string',
					help='determine the SiO2 model: a-Quartz, cristobalite')
	parser.add_option('-p','--periodic',dest='pbc',action='store_true',\
					default=False,help='determine the periodicity of structure')
	(options, args) = parser.parse_args(argv[1:])
    
    #manage parse errors
	if options.credits: #display credits and quit
		credits="\n**********************************\n\
		Kim QHwan\n\
		Contacts: kim.qhwan@gmail.com\
		\n*********************************\n"
		print(credits)
		exit(0)
    
	if len(args)==0:   #arguments missing
		parser.exit(parser.print_help())
    
	if len(args)>1: #check if more than one argument (NOT OPTION) has been parsed
		parser.error('You have given me more than one argument '+str(args)+'... dunno what to do...\n')
    
	if options.model != 'a-quartz' and options.model != 'cristobalite': 
		parser.error('Unknown model: valid structures is aQuartz')
    
	return options, args

def appendLine(natom, coord, basis, nresidue, atomname):
	natom += 1
	coord1 = coord+basis
	tmpcoords=[]
	tmpcoords.append(nresidue)
	tmpcoords.append('CRI')
	tmpcoords.append(str(atomname))
	tmpcoords.append(natom)
	tmpcoords.append(coord1[0])
	tmpcoords.append(coord1[1])
	tmpcoords.append(coord1[2])
	return natom, tmpcoords


def cristobalite(x, y,z):
	nsi = 0
	no = 0
	noh = 0
	nhh = 0
	atc=[]
	nresidue=0
	natom=0
	SI1 = array([0., 0., 0.])
	SH1 = array([-2.46, 1.42, -1.])
	SI7 = array([-2.46, 1.42, -1.])
	SI2 = array([-2.46, 4.26, 0.])
	SH2 = array([0., 5.68, -1])
	SI8 = array([0., 5.68, -1])
	O1 = array([-1.23, 0.71, -0.5])
	O2 = array([0., 0., 1.51])
	O3 = array([1.23, 0.71, -0.5])
	O4 = array([-2.46, 2.84, -0.5])
	O5 = array([1.23, 4.97, -0.5])
	O6 = array([-2.46, 4.26, 1.51])
	O7 = array([-1.23, 4.97, -0.5])
	O8 = array([0, 7.1, -0.5])
	SI3 = array([0., 0., 3.02])
	SI4 = array([-2.46, 4.26, 3.02])
	SI5 = array([-2.46, 1.42, 4.02])
	SI6 = array([0., 5.68, 4.02])
	SH3 = array([-2.46, 1.42, 4.02])
	SH4 = array([0., 5.68, 4.02])
	O9 = array([-1.23, 0.71, 3.52])
	O10 = array([1.23, 0.71, 3.52])
	O11 = array([-2.46, 2.84, 3.52])
	O12 = array([1.23, 4.97, 3.52])
	O13 = array([-1.23, 4.97, 3.52])
	O14 = array([0., 7.1, 3.52])
	OH1 = array([-2.46, 1.42, -2.51])
	OH2 = array([-0, 5.67, -2.51])
	O15 = array([-2.46, 1.42, -2.51])
	O16 = array([-0, 5.67, -2.51])
	OH3 = array([-2.46, 1.42, 5.53])
	OH4 = array([0., 5.68, 5.53])

	# number of unit cell in a row
	x1 = 4.92
	y1 = 8.52
	z1 = 8.04
	a_num=int(x/x1)
	b_num=int(y/y1)
	c_num=int(z/z1)
	a = array([4.92, 0., 0.])
	b = array([0., 8.52, 0.])
	c = array([0., 0., 8.04])
	# index of water molecule in a row
	a_i=-1
	b_i=-1
	c_i=-1
	coord = array([0., 0., 0.])
    #build row for X
	while c_i<c_num:
		c_i+=1
		b_i=-1
		a_i=-1
		while b_i<b_num:
			b_i+=1
			#coord=a2*a2_i+a3*a3_i
			a_i=-1
			while a_i<a_num:
				a_i+=1
				nresidue+=1
				coord = a_i*a + b_i*b + c_i*c
				theta = uniform(0,360)
				dx = cos(theta*3.14/180)
				dy = sin(theta*3.14/180)
				HH1 = array([-2.46+dx, 1.42+dy, -2.84])
				HH2 = array([0.+dx, 5.68+dy, -2.84])
				HH3 = array([-2.46+dx, 1.42+dy, 5.86])
				HH4 = array([0.+dx, 5.68+dy, 5.86])
				# Append the Atoms
				natom, tmpcoords = appendLine(natom, coord, SI1, nresidue, 'SI1')
				atc.append(tmpcoords)
				nsi += 1
				natom, tmpcoords = appendLine(natom, coord, SI2, nresidue, 'SI2')
				atc.append(tmpcoords)
				nsi += 1
				natom, tmpcoords = appendLine(natom, coord, SI3, nresidue, 'SI3')
				atc.append(tmpcoords)
				nsi += 1
				natom, tmpcoords = appendLine(natom, coord, SI4, nresidue, 'SI4')
				atc.append(tmpcoords)
				nsi +=1 
				natom, tmpcoords = appendLine(natom, coord, O1, nresidue, 'O1')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O2, nresidue, 'O2')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O3, nresidue, 'O3')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O4, nresidue, 'O4')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O5, nresidue, 'O5')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O6, nresidue, 'O6')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O7, nresidue, 'O7')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O8, nresidue, 'O8')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O9, nresidue, 'O9')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O10, nresidue, 'O10')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O11, nresidue, 'O11')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O12, nresidue, 'O12')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O13, nresidue, 'O13')
				atc.append(tmpcoords)
				no += 1
				natom, tmpcoords = appendLine(natom, coord, O14, nresidue, 'O14')
				atc.append(tmpcoords)
				no += 1
				if c_i == 0:
					natom, tmpcoords = appendLine(natom, coord, SH1, nresidue, 'SH1')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, SH2, nresidue, 'SH2')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, OH1, nresidue, 'OH1')
					atc.append(tmpcoords)			
					noh += 1
					natom, tmpcoords = appendLine(natom, coord, OH2, nresidue, 'OH2')
					atc.append(tmpcoords)	
					noh += 1
					natom, tmpcoords = appendLine(natom, coord, HH1, nresidue, 'HH1')
					atc.append(tmpcoords)	
					nhh += 1
					natom, tmpcoords = appendLine(natom, coord, HH2, nresidue, 'HH2')
					atc.append(tmpcoords)
					nhh += 1
					natom, tmpcoords = appendLine(natom, coord, SI5, nresidue, 'SI5')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, SI6, nresidue, 'SI6')
					atc.append(tmpcoords)
					nsi += 1
				elif c_i == c_num:
					natom, tmpcoords = appendLine(natom, coord, SH3, nresidue, 'SH3')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, SH4, nresidue, 'SH4')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, OH3, nresidue, 'OH3')
					atc.append(tmpcoords)			
					noh += 1
					natom, tmpcoords = appendLine(natom, coord, OH4, nresidue, 'OH4')
					atc.append(tmpcoords)	
					noh += 1
					natom, tmpcoords = appendLine(natom, coord, HH3, nresidue, 'HH3')
					atc.append(tmpcoords)	
					nhh += 1
					natom, tmpcoords = appendLine(natom, coord, HH4, nresidue, 'HH4')
					atc.append(tmpcoords)
					nhh += 1
					natom, tmpcoords = appendLine(natom, coord, SI7, nresidue, 'SI7')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, SI8, nresidue, 'SI8')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, O15, nresidue, 'O15')
					atc.append(tmpcoords)			
					no += 1
					natom, tmpcoords = appendLine(natom, coord, O16, nresidue, 'O16')
					atc.append(tmpcoords)
					no += 1
				else:
					natom, tmpcoords = appendLine(natom, coord, SI7, nresidue, 'SI7')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, SI8, nresidue, 'SI8')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, O15, nresidue, 'O15')
					atc.append(tmpcoords)			
					no += 1
					natom, tmpcoords = appendLine(natom, coord, O16, nresidue, 'O16')
					atc.append(tmpcoords)
					no += 1
					natom, tmpcoords = appendLine(natom, coord, SI5, nresidue, 'SI5')
					atc.append(tmpcoords)
					nsi += 1
					natom, tmpcoords = appendLine(natom, coord, SI6, nresidue, 'SI6')
					atc.append(tmpcoords)
					nsi += 1
	print('\n***********************************')
	print('SiO2 surface: a=')
	a_pbc=4.9134
	b_pbc=8.51
	print('Periodic (if apply) (ang): a= ',a_pbc, ' b= ',b_pbc)
	print(nsi, no, noh, nhh)
	return atc, a_pbc*(a_num+1), b_pbc*(b_num+1)



def aQuartz(x, y, z):
	atc=[]
	nresidue=1	# this number will be never changed
	natom=0
	# primitive vector
	a=array([4.9134,0.,0.])
	b=array([0.,8.51,0.])
	c=array([0.,0.,5.4052])
	# basis vector -- 1
	a1=array([4.9134,0.,0.])
	a2=array([-2.4567,4.2551,0.])
	a3=array([0.,0.,5.4052])
	# fractional coefficient
	x1=0.4697;	x2=0.4133;	y2=0.2672;	z2=0.1188
	# basis vector -- 2
	b1 = a1*x1 + a3*0.667
	b2 = a2*x1 + a3*0.333
	b3 = a1*(-x1) + a2*(-x1)
	b4 = a1*x2 + a2*y2 + a3*z2
	b5 = a1*(-y2) + a2*(x2-y2) + a3*(0.667+z2)
	b6 = a1*(y2-x2) + a2*(-x2) + a3*(0.333+z2)
	b7 = a1*y2 + a2*x2 - a3*z2
	b8 = a1*(-x2) + a2*(y2-x2) + a3*(0.667-z2)
	b9 = a1*(x2-y2) - a2*y2 + a3*(0.333-z2)
	b10 = b3+a1+a2
	b11 = b2-a2
	b12 = b1-a1-a2
	b13 = b8+a1+a2
	b14 = b9+a2
	b15 = b6+a2
	b16 = b5-a2-a3 #a3 is used for Si-terminated surface
	b17 = b4-a1-a2
	b18 = b7-a1-a2
	# number of unit cell in a row
	x1 = 4.9134
	y1 = 8.5100
	z1 = 5.4052
	a_num=int(x/x1)
	b_num=int(y/y1)
	c_num=int(z/z1)
	# index of water molecule in a row
	a_i=-1
	b_i=-1
	c_i=-1
    #build row for X
	while c_i<c_num:
		c_i+=1
		b_i=-1
		while b_i<b_num:
			b_i+=1
			#coord=a2*a2_i+a3*a3_i
			a_i=-1
			while a_i<a_num:
				a_i+=1
				coord=a*a_i+b*b_i+c*c_i
				xcoord = coord[0]	# these two are used for pbc
				ycoord = coord[1]
				#1. Append the Si1
				natom += 1
				coord1=coord+b1
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('SI1')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
                #2. Append the Si2
				natom += 1
				coord1=coord+b2
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('SI2')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
                #3. Append the Si3
				natom += 1
				coord1=coord+b3
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('SI3')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
                #4. Append the O1
				natom += 1
				coord1=coord+b4
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O1')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
                #5. Append the O2
				natom += 1
				coord1=coord+b5
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O2')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
                #6. Append the O3
				natom += 1
				coord1=coord+b6
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O3')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
                #7. Append the O4
				natom += 1
				coord1=coord+b7
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O4')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
                #8. Append the O1
				natom += 1
				coord1=coord+b8
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O5')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
                #9. Append the O6
				natom += 1
				coord1=coord+b9
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O6')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
				
				#10. Append the SI4
				natom += 1
				coord1=coord+b10
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('SI4')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
				#11. Append the SI5
				natom += 1
				coord1=coord+b11
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('SI5')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
				#12. Append the SI6
				natom += 1
				coord1=coord+b12
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('SI6')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
				#13. Append the O7
				natom+=1
				coord1=coord+b13
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O7')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
				#14. Append the O8
				natom+=1
				coord1=coord+b14
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O8')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
				#15. Append the O9
				natom+=1
				coord1=coord+b15
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O9')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
				#16. Append the O10
				natom+=1
				coord1=coord+b16
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O10')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
				#17. Append the O11
				natom+=1
				coord1=coord+b17
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O11')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
				#18. Append the O12
				natom+=1
				coord1=coord+b18
				tmpcoords=[]
				tmpcoords.append(nresidue)
				tmpcoords.append('AQU')
				tmpcoords.append('O12')
				tmpcoords.append(natom)
				tmpcoords.append(coord1[0])
				tmpcoords.append(coord1[1])
				tmpcoords.append(coord1[2])
				atc.append(tmpcoords)
                    
	print('\n***********************************')
	print('SiO2 surface: a= ',abs(xcoord),' b= ',abs(ycoord))
	a_pbc=4.9134
	b_pbc=8.51
	print('Periodic (if apply) (ang): a= ',a_pbc, ' b= ',b_pbc)
	return atc, a_pbc, b_pbc

def write_gro(file, data, pbc1="",pbc2=""):
	'''
	Write a gromacs gro file.
		Input Variables:
			file: output file (type: file)
			data: list of lists. Each list contains:
				  1) residue number
				  2) residue name(WATER)
				  3) atom name
				  4) atom number
				  5,6,7) X-, Y-, and Z- coordinates
		Variables:
			line: store each list of data(type: list)
			outline: string containing a single line to be witten in file
					 (type: string)
				
	'''
	file.write("Generated by buildSiO2 v1.1\n" + str(len(data)) + "\n")
	for index, line in enumerate(data):
		outline="%5i%5s%5s%5i%8.3f%8.3f%8.3f" % (line[0],line[1],line[2],line[3],float(line[4])/10.0,float(line[5])/10.0,float(line[6])/10.0)
		file.write(outline+"\n")
	outline="	10	10	10\n"
	if pbc1 =="":
		outline="	10	10	10\n"
	elif pbc2 =="":
		outline="	10	"+str(float(pbc1)/10.0)+"	10\n"
	else:
		outline="	"+str(float(pbc1)/10.0)+"	"+str(float(pbc2)/10.0)+"	10\n"
	file.write(outline+"\n")
	return


def main():
#===================================================#
#				MAIN	MAIN	MAIN				#
#===================================================#
	(options,args)=parsecmd()

	ofile=args[0]	#get output file to save structure

	if options.model == 'a-Quartz': 
		atc, pbc_a, pbc_b=aQuartz(options.box[0],options.box[1],options.box[2])
	elif options.model == 'cristobalite':
		atc, pbc_a, pbc_b=cristobalite(options.box[0], options.box[1], options.box[2])
	print('saving structure...')
	print('******************************')
	OUT=open(ofile,'w')
	if options.pbc:
		write_gro(OUT,atc,pbc_a,pbc_b)
	else:
		write_gro(OUT,atc)
	exit(0)

if __name__ == "__main__":
	main()
