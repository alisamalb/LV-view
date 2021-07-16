import numpy as np
class Structure:
	def __init__(self):
		self.atomlines=0
		self.atoms=[]
		self.coord=''
	def addAtom(self,resname,name,x,y,z):
		if len(self.atoms)>0:
			resnr=self.atoms[-1].resnr+1
		else:
			resnr=1
		number=(len(self.atoms)+1)%99999
		self.atoms.append(Atom("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (resnr,resname,name,number,x,y,z)))
	def write(self, path):
		f=open(path,"w")
		f.write("GROMACS\n")
		f.write(str(len(self.atoms))+"\n")
		for a in self.atoms:
			f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (a.resnr,a.resname,a.name,a.number,a.x,a.y,a.z) )
		f.write(self.coord)
		f.close()
			
class Atom:
	def __init__(self,line):
		self.resnr=int(line[:5])%99999
		self.resname=line[5:10]
		self.name=line[10:15]
		self.number=int(line[15:20])
		self.x=float(line[20:28])
		self.y=float(line[28:36])
		self.z=float(line[36:44])

def fixed_string(string,length):
	return " "*(length-len(str(string)))+str(string)
def read_index(path):
	'''Imports a .ndx files and returns a 2-row array of the group names and atom IDs'''
	file=open(path,"r").readlines()
	file=''.join(' '.join(file).split('\n'))
	file=file.split('[')
	file=[x for x in file if len(x)>0]
	names=[x.split(']')[0] for x in file]
	index=[x.split(']')[1] for x in file]
	index=[x.split(' ') for x in index]
	index=[[int(y) for y in x if len(y)>0] for x in index]
	return [names,index]

def read_structure(path):
	'''reads a .gro and outputs a Structure object'''
	file=open(path,"r").readlines()
	file=[x for x in file if len(x)>0] #removes empty lines
	structure=Structure()
	structure.atoms=[Atom(x) for x in file[2:-1]]
	structure.atomlines=int(file[1])
	structure.coord=file[-1]
	return structure

def geom_center(atoms, index):
	'''calculates the geometric center of the atoms selectted from the "atoms" list using the index'''
	sel_atoms=[atoms[x-1] for x in index]
	sum_x=sum([a.x for a in sel_atoms])/len(sel_atoms)
	sum_y=sum([a.y for a in sel_atoms])/len(sel_atoms)
	sum_z=sum([a.z for a in sel_atoms])/len(sel_atoms)
	return (sum_x, sum_y, sum_z)

def calculate_axis(atoms,index):
	print("-- LV-VIEW --")
	print("Select molecule A")
	for i in range(len(index[0])):
		print(str(i)+". "+index[0][i])
	mol_A=int(input("\nSelect group: "))
	print("----\nSelect molecule B")
	for i in range(len(index[0])):
                print(str(i)+". "+index[0][i])
	mol_B=int(input("\nSelect group: "))
	A_center=geom_center(atoms,index[1][mol_A])
	B_center=geom_center(atoms,index[1][mol_B])
	return (A_center,B_center)

def write_sphere(structure, rho, maxtau,COM1):
	j=np.linspace(-rho,rho,int(rho*50))
	for x in j:
		for y in j:
			if not (np.sqrt(x*x+y*y)>rho):
				for sign in [-1,1]:
					z=sign*np.sqrt(-x*x-y*y+rho*rho)
					position=(x+COM1[0],y+COM1[1],z+COM1[2])
					theta=np.arctan2(z,y)
					distance=np.sqrt(x*x+y*y+z*z)
					tau=np.sqrt(x+distance)
					taulim=maxtau
					if  tau<taulim:
						structure.addAtom("DUM","DUM",position[0],position[1],position[2])
					if taulim-0.005<tau<taulim+.005:
						for i in range(int(rho)*10):
							structure.addAtom("DUM","DUM",x*(i/rho/10)+COM1[0],y*np.sqrt(i/rho/10)+COM1[1],np.sqrt(i/rho/10)*z+COM1[2])

def suggest_rotation(COM1,COM2):
	vector=[COM1[0]-COM2[0],COM1[1]-COM2[1],COM1[2]-COM2[2]]		
	print("The COM1 -> COM2 vector is {}".format(vector))
	dist=np.sqrt(vector[0]**2+vector[1]**2)
	rot=-np.arcsin(vector[1]/dist)*180/3.14
	print("The suggested rotation is 0 0 "+str(rot)+" degrees.")
	print("The initial tau angle is "+str(np.sqrt(vector[0]+dist)))
index=read_index(input("Path to the index file: "))
structure=read_structure(input("Path to the structure: "))
COM1,COM2=calculate_axis(structure.atoms,index)
suggest_rotation(COM1,COM2)
print("The COM1 is {}".format(COM1))
rho=float(input("Max radius of the sphere (rho): "))
maxtau=float(input("Max tau angle: "))
volume=Structure()
volume.coord=structure.coord
write_sphere(volume,rho,maxtau,COM1)
volume.write("volume.gro")

