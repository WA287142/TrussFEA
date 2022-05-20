from dis import dis
from errno import EACCES
from json import load
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
numNodes = 0
numElems = 0
nodes = {}
elemDict = {}
ymod = 0
area = 0
EA = 0
numCons = 0
NOTdofs = []
pForces = 0 #number of P forces
pForceDict = {}
qForces = 0
qForceDict = {}
nodeDOFs = {}
globalstiff = np.zeros([2*numNodes, 2*numNodes])
loadVector = []
nodaldisplacements = {}
eleminfo = []
nodeinfo = []
q_sol = []
def readInputs(fname, nodes):
f = open(fname)
global numNodes
global numElems
global eleminfo
line = f.readline().split()
numNodes = int(line[0])
numElems = int(line[1])
for i in range(numNodes): # reads the lines with node #s on it and
puts them into a node and xcoords array
line = f.readline().split()
nodes[int(line[0])] = [float(line[1]), float(line[2])]
print("Nodes and X,Y Coords: ")
for key in nodes:
print(key, " : ", nodes[key])
coords = nodes[key]
temp = [key, int(coords[0]), int(coords[1])]
nodeinfo.append(temp)
print("Node info: ", nodeinfo)
f.readline() # gets rid of the blank line
# Reading the elements
for i in range(numElems):
line = f.readline().split()
elemDict[int(line[0])] = [line[1], line[2]]
print("\nelement dictionary: ")
for key in elemDict:
print(key, " : ", elemDict[key])
connects = elemDict[key]
temp = [key, int(connects[0]), int(connects[1])]
eleminfo.append(temp)
print("ELem info: ", eleminfo)
print("number of elements: ", numElems)
f.readline() # Gets rid of the blank line
# Reading in Youngs modulus and area A
line = f.readline().split()
global ymod
temp = line[0].split('e')
ymod = float(temp[0]) * 10**float(temp[1])
global area
area = float(line[1])
global EA
EA = ymod * area
f.readline() # Gets rid of blank line
#Number of Constraints
global numCons
global NOTdofs
numCons = int(f.readline())
NOTdofs = f.readline().split()
# print("Number of Constraints: ", numCons, " \nList of constrained
DOFs are: ", NOTdofs)
# Reading in the Horizontal and Vertical Forces
global pForces
f.readline() #get rid of blank line
line = f.readline() #P forces are vertical
pForces = int(line)
# print("# of P forces: ", pForces)
for i in range(pForces):
line = f.readline().split()
pForceDict[int(line[0])] = float(line[1])
# print("P forces: ",pForceDict)
global qForces
f.readline() #get rid of blank line
line = f.readline()
# print("line:", line)
qForces = int(line)
# print("# of Q forces: ", qForces)
for i in range(qForces):
line = f.readline().split()
qForceDict[int(line[0])] = float(line[1])
def getElementK(element):
# localstiff = np.zeros([4,4])
elemNodes = elemDict[element]
node1 = nodes[int(elemNodes[0])]
node2 = nodes[int(elemNodes[1])]
node1x = node1[0]
node1y = node1[1]
node2x = node2[0]
node2y = node2[1]
l = np.sqrt((node2y - node1y)**2 + (node2x - node1x)**2)
theta = np.arcsin((node2y - node1y)/l)
#print("Theta (rads):", theta)
#print("Length: ", l)
# K[0,0] = 1 * EA / l
# K[0,1] = -1 * EA / l
# K[1,0] = -1 * EA / l
# K[1,1] = 1 * EA / l
eal = EA / l
trussK = np.zeros([4,4])
trussK[0,0] = eal * np.cos(theta)**2
trussK[0,1] = eal * np.cos(theta) * np.sin(theta)
trussK[0,2] = -eal * np.cos(theta)**2
trussK[0,3] = -eal * np.cos(theta) * np.sin(theta)
trussK[1,0] = eal * np.cos(theta) * np.sin(theta)
trussK[1,1] = eal * np.sin(theta)**2
trussK[1,2] = -eal * np.cos(theta) * np.sin(theta)
trussK[1,3] = -eal * np.sin(theta)**2
trussK[2,0] = -eal * np.cos(theta)**2
trussK[2,1] = -eal * np.cos(theta) * np.sin(theta)
trussK[2,2] = eal * np.cos(theta)**2
trussK[2,3] = eal * np.cos(theta) * np.sin(theta)
trussK[3,0] = -eal * np.cos(theta) * np.sin(theta)
trussK[3,1] = -eal * np.sin(theta)**2
trussK[3,2] = eal * np.cos(theta) * np.sin(theta)
trussK[3,3] = eal * np.sin(theta)**2
# print("Truss K = \n", trussK)
return trussK
def getElementDOF(element):
global nodeDOFs
# print("\nGetting DOFs for element", element, "\n")
# input parameter element should be the element #, an int
# two dofs per NODE, not element
# Find all dofs per node first
for i in nodes:
dof1 = int(i) * 2 - 1
dof2 = int(i) * 2
nodeDOFs[int(i)] = dof1, dof2
# print("dofs for node", i, ": ", dof1, dof2)
# print("nodes for element", element, ":", elemDict[element])
# print("DOFs that are 0:", dofs)
# print("element", element, "DOFs:", elemDOFs)
elemDOFs = []
nodeValues = elemDict[element] # stores the node values in
nodeValues
temp = nodeDOFs[int(nodeValues[0])] # Stores the dofs of the first
node in temp
elemDOFs.append(int(temp[0]))
elemDOFs.append(int(temp[1]))
temp = nodeDOFs[int(nodeValues[1])] # Stores the dofs of the second
node in temp
elemDOFs.append(int(temp[0]))
elemDOFs.append(int(temp[1]))
# print("nodes for element", element, ":", elemDict[element])
# print("Element", element, "DOF values: ", elemDOFs)
return elemDOFs
def assembleGlobalStiffnessMatrix():
global globalstiff
globalstiff = np.zeros([2*numNodes, 2*numNodes])
for i in range(numElems):
elemDOFs = getElementDOF(i+1)
# print("ElemDOFs for element", i, ": ", elemDOFs)
elemK = getElementK(i+1)
for row in range(len(elemDOFs)):
globalRow = elemDOFs[row] - 1
for col in range(len(elemDOFs)):
globalCol = elemDOFs[col] - 1
#print("stifness for:", globalRow, globalCol)
globalstiff[globalRow, globalCol] = globalstiff[globalRow,
globalCol] + elemK[row,col]
print("\nGlobal Stiffness Matrix after element", i+1, " \n")
print(globalstiff)
# Load Vector
global loadVector
loadVector = np.zeros([2*numNodes,1])
for i in pForceDict:
dof = (i * 2)
loadVector[dof - 1] += pForceDict[i]
for i in qForceDict:
dof = (i * 2) - 1
loadVector[dof - 1] += qForceDict[i]
print("Load Vector:\n", loadVector)
def imposeConstraints():
tempstiff = globalstiff
for i in NOTdofs:
globalstiff[int(i)-1] = 0
globalstiff[:, int(i)-1] = 0
globalstiff[int(i)-1, int(i)-1] = 1
# globalstiff[int(i)-1, int(i)-1] = globalstiff[int(i)-1,
int(i)-1] + (1 * 10**30)
print("DOFs that need to be constrained:", NOTdofs)
print("Global Stiffness Matrix after imposed constraints: \n",
globalstiff)
def solver():
global nodaldisplacements
q = np.linalg.solve(globalstiff, loadVector)
# print(q)
j = 1
for i in range(0, 2 * numNodes - 1, 2):
nodaldisplacements[j] = [float(q[i]), float(q[i+1])]
j += 1
# Applying the moment
# momentVec = np.zeros([2*numNodes,1]) # all the DOFs
# pointForce[momentDOF-1] = moment
print("\nNodal Displacements:")
for node in nodaldisplacements:
print(node, " : ", nodaldisplacements[node])
temp = nodaldisplacements[node]
q_sol.append(temp[0])
q_sol.append(temp[1])
# print("q_sol = ", q_sol)
return nodaldisplacements
# print("nodal displacements:\n", nodaldisplacements)
def calulateLocalForce(element):
# To start, apply the nodal displacements to the node coordinates
# Put into a new node dict
elemNodes = elemDict[element]
node1 = nodes[int(elemNodes[0])] # node1 is a list of x y
coords, elemNodes[0] is the first node
node2 = nodes[int(elemNodes[1])]
node1x = node1[0]
node1y = node1[1]
node2x = node2[0]
node2y = node2[1]
l = np.sqrt((node2y - node1y)**2 + (node2x - node1x)**2)
theta = np.arcsin((node2y - node1y)/l) # theta of
element
tRot = np.zeros([4,4])
tRot[0,0] = np.cos(theta)
tRot[0,1] = np.sin(theta)
tRot[1,0] = -np.sin(theta)
tRot[1,1] = np.cos(theta)
tRot[2,2] = np.cos(theta)
tRot[2,3] = np.sin(theta)
tRot[3,2] = -np.sin(theta)
tRot[3,3] = np.cos(theta)
# print("Theta: ", theta)
# print("T rotation: \n", tRot)
tDim = np.zeros([2, 4])
tDim[0,0] = 1 #T dimensional, [[1 0 0 0], [0
0 1 0]]
tDim[1,2] = 1
# print("TDim = ", TDim)
trussT = tDim @ tRot # tRot is a 4x4, TDim is 2x4 =
2x4
q1 = nodaldisplacements[int(elemNodes[0])] # elemNodes[0] is the
1st node, disp1 is the deflection at that node
q_x1 = q1[0]
q_y1 = q1[1]
q2 = nodaldisplacements[int(elemNodes[1])] # elemNodes[1] is the
2nd node, disp1 is the deflection at that node
q_x2 = q2[0]
q_y2 = q2[1]
qTruss = np.zeros([4,1]) # global hori and vert
deflections of element e , 4x1 matrix
qTruss[0,0] = q_x1
qTruss[1,0] = q_y1
qTruss[2,0] = q_x2
qTruss[3,0] = q_y2
# print("qTruss: \n", qTruss)
q_e = trussT @ qTruss # trussT 2x4 and qTruss a
4x1 makes a 2x1
# print("q_e = \n", q_e)
# Q_e should be 1x4 or 4x1, K_e @ q_e -> 1x4, q_e is 2x1 so K_e has to
be a 4x2 but is a 4x4 right now
u1 = q_e[0]
u2 = q_e[1]
ax_strain = -(u1 * 1/l) + (u2 * 1/l)
# print("axial strain: ", ax_strain)
ax_stress = ax_strain * ymod
print("Element stress at two ends for element ", element, " : ",
ax_stress, ax_stress)
# force = ax_stress * ax_strain
# print("force = ", force)
return ax_stress, ax_stress
def show_truss(q_sol, elem_num, elems_info, nodes_info):
''' function entries:
1. q_sol: LIST of solutions to global nodal displacements
- format: [qx_node1, qy_node1, qx_node2, qy_node2,
...,qx_nodeN, qy_nodeN],
where; - qx_nodeN represents displacement of node N in
x-direction
- qy_nodeN represents displacement of node N in
y-direction
2. elem_num: number of elements in truss
- format: An integer number
3. elems_info: element connectivity information
- format: [[elem1, elem1_node1, elem1_node2], [elem2,
elem2_node1, elem2_node2], ...,[elemN, elemN_node1, elemN_node2]]
where; - elemN represents the element number N (an
integer)
- elemN_node1 represents node1 of element N (an
integer)
- elemN_node2 represents node2 of element N (an
integer)
4. nodes_info: node global coordinates
- format: [[node1, x1_coord, y1_coord], [node2, x2_coord,
y2_coord], ...,[nodeN, xN_coord, yN_coord]]
where; - nodeN represents node number N (an integer)
- xN_coord represents x location of node N
- yN_coord represents y location of node N
For Example: JDW Example - Chapter 4, page 14;
elems_info = [[1, 1, 2], [2, 1, 4], [3, 3, 2], [4, 4, 2],
[5, 3, 4]]
nodes_info = [[1, 0, 10], [2, 10, 10], [3, 0, 0], [4, 10,
0]]
q_sol = [0, 0, -0.250000000000000, 0.957106781186548, 0, 0,
0, 0.707106781186548]'''
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#-------------------------------------------------------------------------
figt, tplt = plt.subplots(1, 1)
for e in range(elem_num):
e_info = elems_info[e]
e_nodes = [e_info[1], e_info[2]]
#---------------------------------------------------------------------
# getting element dofs
e_dofs = []
for i in range(len(e_nodes)):
e_dofs.append(2*e_nodes[i] - 1)
e_dofs.append(2*e_nodes[i])
#---------------------------------------------------------------------
# computing node length
node1_info, node2_info = nodes_info[e_nodes[0]-1],
nodes_info[e_nodes[1]-1]
node1_x, node1_y, node2_x, node2_y = node1_info[1], node1_info[2],
node2_info[1], node2_info[2]
le_final = ((node2_x - node1_x)**2 + (node2_y - node1_y)**2)**0.5
le_final = float(le_final)
#---------------------------------------------------------------------
# getting element nodal displacements
qe_final = []
for dof in e_dofs:
qe_final.append(q_sol[dof - 1])
#print(dof, q_sol[dof - 1])
#---------------------------------------------------------------------
# deformed_coord = initial_coord + deformation, qe
x_undef_range, y_undef_range = [node1_x, node2_x], [node1_y,
node2_y]
x_def_range, y_def_range = [node1_x+float(qe_final[0]),
node2_x+float(qe_final[2])], [node1_y+float(qe_final[1]),
node2_y+float(qe_final[3])]
#---------------------------------------------------------------------
# Plots
if e != elem_num-1:
tplt.plot(x_undef_range, y_undef_range, 'g--')
tplt.plot(x_def_range, y_def_range, color='black')
else:
tplt.plot(x_undef_range, y_undef_range, 'g--',
label='Undeformed')
tplt.plot(x_def_range, y_def_range, color='black',
label='Deformed')
tplt.set(xlabel='x, mm', ylabel='y, mm')
tplt.set_title('Truss Deformation')
tplt.grid()
tplt.legend()
plt.show()
return figt
def ReportResults(fname):
global numNodes
global nodes
global numElems
global elemDict
global ymod
global area
global numCons
global load
global NOTdofs
global pForceDict
global qForceDict
f = open(fname + "output.txt", 'w')
f.write("Number of Nodes: \n")
f.write(str(numNodes))
f.write("\nNumber of Elements: \n")
f.write(str(numElems))
f.write("\nNodes and X,Y Coords: \n")
for key in nodes:
line = str(key) + " : " + str(nodes[key])
f.write(line)
f.write("\nElement dictionary: \n")
for key in elemDict:
line = str(key) + " : " + str(elemDict[key])
f.write(line)
f.write("\nYoung's Modulus/ EI: \n")
f.write(str(ymod))
f.write("\nNumber of Constraints: \n")
f.write(str(numCons))
f.write("\nList of DOFs that are constrained: \n")
f.write(str(NOTdofs))
f.write("\nVertical forces: " + str(pForceDict))
f.write("\nHorizontal forces: " + str(qForceDict))
for i in range(numElems):
elemDOFs = getElementDOF(i+1)
line = "\nElement DOFs for element " + str(i) + ": \n" +
str(elemDOFs)
f.write(line)
elemK = getElementK(i+1)
line = "\nLocal Stiffness Matrix for element " + str(i) + ": \n" +
str(elemK)
f.write(line)
f.write("\nGlobal Stiffness Matrix: \n")
assembleGlobalStiffnessMatrix()
f.write(str(globalstiff))
imposeConstraints()
f.write("\nGlobal Stiffness Matrix after imposing constraints: \n")
f.write(str(globalstiff))
f.write("\nNodal Displacements:\n")
for node in nodaldisplacements:
f.write(str(node) + " : " + str(nodaldisplacements[node]) + "\n")
for i in range(numElems):
line = calulateLocalForce(i+1)
f.write("\nAxial stress for element " + str(i) + " = " +
str(line))
readInputs("Project2.Bonus.txt", nodes)
assembleGlobalStiffnessMatrix()
imposeConstraints()
solver()
show_truss(q_sol, numElems, eleminfo, nodeinfo)
for i in range(numElems):
calulateLocalForce(i+1)
ReportResults("Project2.Bonus.txt")