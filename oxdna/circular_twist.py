from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


def get_array(line):   #turns a line in a trajectpory file into a numpy array
    return np.array([float(line[0][:-1]), float(line[1][:-1]), float(line[2])])

def vector_string(array):
    return str(array[0])+", "+str(array[1])+", "+str(array[2])



#read and write appropriate files
init = open("circ_300.dat", "r").readlines()
input_file = open("input.dat", "w")
external_forces = open("forces.dat","w")


N_traps = 2 #number of PAIRS of twist traps to be generated
rate = 1e-5
stiff = 11
rate = 1e-5
base = 0.
group_name = ""
mask = [1,1,1]


#turns lists which are read from file into numpy arrays
def vectorize(thing):
    return np.array([float(thing[0]),float(thing[1]),float(thing[2])])


N_particles = int(len(init)-3) #total number of particles
N = int(N_particles/2)              #number of particles in each signle strand
ratio = int(N/N_traps)                    #ratio of number of praticles to each trap

def get_centre(A, B):       #find the centre point of two quantities
    return 0.5*(A+B)


list_of_centres=[]
positions=[]
compl_positions=[]

for i in range(N):
    traj_line = get_array(init[3+i].split())      #particle i position is in line i+3
    comp_traj_line = get_array(init[3+2*N-1-i].split())
    particle = (traj_line)
    complement  = comp_traj_line      #find position of other particle in bp
    list_of_centres.append(get_centre(particle,complement))     #finds centre of these two points and adds them to list
    positions.append(particle)  
    compl_positions.append(complement)        


circle_centre = get_centre(list_of_centres[0], list_of_centres[N-1])   #find the centre of the circle

circle_normal = np.cross(list_of_centres[0]-circle_centre, list_of_centres[1]-circle_centre)  #find the normal to the circular plae
circle_normal = circle_normal/(np.dot(circle_normal,circle_normal))**0.5   


#suppose we want N_traps number of twist traps spread evenly over the circle, thenwe need a axis for each to rotate about
#these axes will lie tangent to the circle and have positions coinciding with the centre of the two nucleotides

tangents = []
for i in range(N_traps):
    tan = np.cross(list_of_centres[ratio*i]-circle_centre, circle_normal)
    tan = tan/(np.dot(tan,tan))**0.5
    tangents.append(tan)
   
#Now we add N_traps numbner of pair of traps    
    
    
  
class ForceTwist:
	def __init__(self, particle, pos0, stiff, rate, base, axis, center, mask, group_name):
		self.particle = str(particle)
		self.pos0 = vector_string(pos0)
		self.stiff = str(stiff)
		self.rate = str(rate)
		self.base = str(base)
		self.axis = vector_string(axis)
		self.center = vector_string(center)
		self.mask = vector_string(mask)
		self.group_name = str(group_name)
	def __str__(self):
		return str(
"{"+"\n"+
"  type = twist\n"+
"  particle = "+self.particle+"\n"+
"  pos0 = "+self.pos0+"\n"+
"  stiff = "+self.stiff+"\n"+ 
"  rate = "+self.rate+"\n"+
"  base = "+self.base+"\n"+
"  axis = "+self.axis+"\n"+
"  center = "+self.center+"\n"+
"  mask = "+self.mask+"\n"+
"  group_name = "+self.group_name+"\n"+
"}"+"\n"+"\n")  
                                                          

#p = ForceTwist(1,[1,1,1], stiff, rate, base , [0,0,0], [0,0,1], mask, group_name)   
#print(p.__str__())
                                                          
twists = []
for i in range(N_traps):
    j=i*ratio
    print(j)
    twists.append(ForceTwist(j, positions[j], stiff, rate, base, tangents[i], list_of_centres[j] , mask, group_name))

  
for i in range(len(twists)):
    external_forces.write(twists[i].__str__())
external_forces.close()







#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
#x_vals = []  
#y_vals = []
#z_vals = []
#for j in range(len(list_of_centres)):
#    x_vals.append(list_of_centres[j][0])
#    y_vals.append(list_of_centres[j][1])
#    z_vals.append(list_of_centres[j][2])
#
#
#
#ax.scatter(x_vals, y_vals, z_vals)
#plt.show()


        

  
 
