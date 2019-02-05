from math import *
import random
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.backend_bases as base
import time
import numpy
import sys


plt.rcParams['toolbar'] = 'None'

def Average(list, n):
    total1 = 0
    total2 =0
    total3 = 0
    for item in list:
        total1 += item[0]
        total2 += item[1]
        total3 += item[2]
    return [total1/n, total2/n, total3/n]

def Roundup(x,a):
  return int(ceil(x / a)) * a

def Norm(vector):
    return sqrt(((abs(vector[0]))**2)+((abs(vector[1]))**2)+((abs(vector[2]))**2))


def timer():
    start = time.time()
    s = System(625, 16, 1, 2.5, 0.01)
    s.initVelocities()
    s.computeForces()
    s.integrateEOM()
    s.run(500)
    end = time.time()
    print (end - start)

def BinCounts(data, width):
    num_bins = int(len(data)/width)
    hist, edges = numpy.histogram(data, bins=num_bins, range=(0,int(width*num_bins)), density=False)
    counts = []
    for i in hist:
        counts.append(i)
    return(counts)

def netForce(forces):
    net = []
    for i in forces:
        net.append(sqrt(i[0]**2+i[1]**2+i[2]**2))
    return net


def potentialLJ(r):
    return 4 * ((r ** -12) - (r ** -6))

exitRequest = False
# some global data that needs to be saved:
fig = plt.figure()
ax1 = plt.gca(projection='3d')
counter = 0
anglecounter = 0
class System:
    def __init__(self, nParticles, boxlength, temperature, rCut, dt, colors, particleSize, bg, cop):
        self.nParticles = nParticles
        self.boxlength = boxlength
        self.prevCoords = []
        self.coords = []
        self.newCoords = []
        self.velocities = []
        self.KE = 0
        self.P = 0
        self.forces = []
        self.temperature = temperature
        self.rCut = rCut
        self.dt = dt
        self.colors = colors
        self.particleSize = particleSize
        self.bg = bg
        self.ani = animation.FuncAnimation(fig, self.animate_, interval=50)
        self.mng = plt.get_current_fig_manager()
        self.cop = cop
        self.dtError = False
        perSide = int(ceil((nParticles)**(1./3.)))
        # Evenly distribute n particles throughout a 3D box with side length=boxlength
        for i in range(1, perSide + 1):
            for j in range(1, perSide + 1):
                for k in range(1, perSide +1):
                    if len(self.coords) < nParticles:
                        self.coords.append(
                            [round((i * ((float(boxlength)) / float(perSide)) - 0.5), 1),
                             round((j * ((float(boxlength)) / float(perSide)) - 0.5), 1),
                             round((k * ((float(boxlength)) / float(perSide)) - 0.5), 1)])

    def initVelocities(self):
        nFreedom = 3*self.nParticles-3

        # Give random velocities to each particle - represented as a list of lists [[x1,y1]...]
        for i in range(self.nParticles):
            self.velocities.append([random.uniform(-0.5, 0.5), random.uniform(-0.5, 0.5),
                                    random.uniform(-0.5, 0.5)])

        # Remove the center-of-mass motion of the system
        sumV = Average(self.velocities, self.nParticles)
        for velocities in self.velocities:
            velocities[0] -= sumV[0]
            velocities[1] -= sumV[1]
            velocities[2] -= sumV[2]

        # Compute kinetic energy and rescale velocity to match temperature
        NormList = []
        for velocities in self.velocities:
            NormList.append(0.5*(Norm(velocities))**2)
        self.KE = sum(NormList)

        scale = sqrt(nFreedom*self.temperature/self.KE)
        for velocities in self.velocities:
            for component in velocities:
                component *= scale

        # Establish the previous position of each particle, according to the time step dt
        for i in range(self.nParticles):
            self.prevCoords.append(
                [self.coords[i][0]-self.velocities[i][0]*self.dt,
                 self.coords[i][1]-self.velocities[i][1]*self.dt,
                 self.coords[i][2] - self.velocities[i][2] * self.dt])


    def computeForces(self):

        # Initialize a force list
        self.forces = []
        for i in range(self.nParticles):
            self.forces.append([0,0,0])

        # Update a list of forces on each particle according to LJ equation ==> indices matches self.coords
        for i in range(self.nParticles):
            for j in range(i+1, self.nParticles):
                distVec = [self.coords[i][0]-self.coords[j][0],
                           self.coords[i][1]-self.coords[j][1],
                           self.coords[i][2] - self.coords[j][2]]
                distVec[0] -= self.boxlength*round(distVec[0]/self.boxlength)
                distVec[1] -= self.boxlength*round(distVec[1]/self.boxlength)
                distVec[2] -= self.boxlength*round(distVec[2]/self.boxlength)
                dist = Norm(distVec)
                if dist == 0:
                    self.dtError = True
                    self.setExit()
                elif dist<=self.rCut:
                    pairwiseForce = (48.*(dist**-14)-(24.*(dist**-8)))
                    self.forces[i][0] += pairwiseForce * distVec[0]
                    self.forces[i][1] += pairwiseForce * distVec[1]
                    self.forces[i][2] += pairwiseForce * distVec[2]
                    self.forces[j][0] -= pairwiseForce * distVec[0]
                    self.forces[j][1] -= pairwiseForce * distVec[1]
                    self.forces[j][2] -= pairwiseForce * distVec[2]
        return(self.forces)

    def integrateEOM(self):

        # Figure out the new coordinates according to Newtonian laws of motion // Verlet algorithm
        for i in range(self.nParticles):
            self.newCoords.append([2.*self.coords[i][0]-self.prevCoords[i][0]+(self.dt**2)*self.forces[i][0],
                             2.*self.coords[i][1]-self.prevCoords[i][1]+(self.dt**2)*self.forces[i][1],
                                  2. * self.coords[i][2] - self.prevCoords[i][2] + (self.dt ** 2) * self.forces[i][2]])
        for element in self.newCoords:
            for component in element:
                component %= self.boxlength

    def update(self):

        # Time steps: what was old is now lost, what was is now old, what was new now is
        self.prevCoords = self.coords
        self.coords = self.newCoords
        self.newCoords = []
        self.computeForces()
        self.integrateEOM()

    def animate_(self, i):
        global counter, anglecounter
        counter += 1

        # Function to called by matlabplot animation (need only have parameter i)
        xs, ys, zs = zip(*self.coords)
        ax1.clear()
        #ax1.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        #ax1.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        #ax1.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax1.w_xaxis.line.set_color((0.0, 0.0, 0.0, 0.0))
        ax1.w_yaxis.line.set_color((0.0, 0.0, 0.0, 0.0))
        ax1.w_zaxis.line.set_color((0.0, 0.0, 0.0, 0.0))
        ax1.grid(False)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_zticks([])
        anglecounter+= 0.75
        #ax1.view_init(15, anglecounter%360)
        if self.cop == 'force':
            ax1.scatter(xs, ys, zs, c=netForce(self.forces), cmap=self.colors, s=self.particleSize)
        elif self.cop == 'velocity':
            ax1.scatter(xs, ys, zs, c=netForce(self.velocities), cmap=self.colors, s=self.particleSize)
        else:
            ax1.scatter(xs,ys,zs, c=self.colors(numpy.linspace(0,1,self.nParticles)), s=self.particleSize)
        ax1.set_facecolor(self.bg)
        self.update()

    def run(self, i):

        # Will just run updates - no animation (useful for timing)
        for i in range(i):
            self.update()

    def animate(self):

        # Function to generate animated simulation
        #plt.subplots_adjust(left=0, bottom=0, right=1, top=1)
        plt.show(block=True)



    def calcPE(self):

        # Calculates PE according to coordinates at any time step
        self.PE = 0
        Ecut = 4*((self.rCut**-12)-(self.rCut**-6))
        for i in range(self.nParticles-1):
            for j in range(i+1, self.nParticles):
                distVec = [self.coords[i][0] - self.coords[j][0],
                           self.coords[i][1] - self.coords[j][1],
                           self.coords[i][2] - self.coords[j][2]]
                distVec[0] -= self.boxlength * round(distVec[0] / self.boxlength)
                distVec[1] -= self.boxlength * round(distVec[1] / self.boxlength)
                distVec[2] -= self.boxlength * round(distVec[2] / self.boxlength)
                dist = Norm(distVec)
                if dist<=self.rCut:
                    self.PE +=4*((dist**-12)-(dist**-6))-Ecut
        return self.PE

    def calcKE(self):

        # Calculates KE according to coordinates at any time step
        dx = []
        self.KE = 0
        for i in range(self.nParticles):
            dx = [self.coords[i][0]-self.prevCoords[i][0], self.coords[i][1]-self.prevCoords[i][1],
                  self.coords[i][2]-self.prevCoords[i][2]]
            dx[0] -= self.boxlength*round(dx[0]/self.boxlength)
            dx[1] -= self.boxlength*round(dx[1]/self.boxlength)
            dx[2] -= self.boxlength*round(dx[2]/self.boxlength)
            self.KE+= 0.5*(Norm(dx)/self.dt)**2
        return self.KE

    def calcP(self):
        # Calculates Pressure according to coordinates and forces
        summation = 0
        density = self.nParticles/(self.boxlength)**3.
        for i in range(self.nParticles-1):
            for j in range(i+1, self.nParticles):
                distVec = [self.coords[i][0] - self.coords[j][0],
                           self.coords[i][1] - self.coords[j][1],
                           self.coords[i][2] - self.coords[j][2]]
                distVec[0] -= self.boxlength * round(distVec[0] / self.boxlength)
                distVec[1] -= self.boxlength * round(distVec[1] / self.boxlength)
                distVec[2] -= self.boxlength * round(distVec[2] / self.boxlength)
                dist = Norm(distVec)
                if dist <=self.rCut:
                    pairwiseForce = (48. * (dist ** -14) - (24. * (dist ** -8)))
                    summation += dist*pairwiseForce
        self.P = float(self.temperature)*density + (1/(3.*self.boxlength**3.))*summation
        return self.P

    def calcGInf(self):

        # Calculates high frequency modulus coefficient GInf re: Zwanzig and Mountain (1965)
        density = self.nParticles / (self.boxlength) ** 3
        return (26/5)*density*self.temperature+3*self.P-(24/5)*density*abs(self.KE)


    def setcoords(self, coords):

        # For experiments: Initialize and then set coordinates to what you want
        for i in range(len(coords)):
            self.coords[i] = coords[i]

    def setExit(self):
        global fig
        global ax1
        self.ani.event_source.stop()
        del self.ani
        plt.close()
        #del self
        fig = plt.figure()
        ax1 = plt.gca(projection='3d')


def runSimulation(s, dt, nRun, nSample, nSampleCoords):
    totalTime = 0
    PElist = []
    KElist = []
    TElist = []
    coordlist = []
    Plist = []
    Glist = []

    for i in range(totalTime, nRun):
        s.computeForces()
        s.integrateEOM()
        s.update()
        if totalTime%nSample == 0:
            PElist.append([totalTime*dt, s.calcPE()])
            KElist.append([totalTime*dt, s.calcKE()])
            TElist.append([totalTime*dt, PElist[-1][1]+KElist[-1][1]])
            Plist.append([totalTime*dt, s.calcP()])
            Glist.append([totalTime*dt, s.calcGInf()])
        if totalTime%nSampleCoords == 0:
            coordlist.append([totalTime*dt, s.coords])
        totalTime+=1
        print(totalTime)
    return PElist, KElist, TElist, coordlist, Plist, Glist

def equilibrate(s, nEquil):
    for i in range(nEquil):
        s.computeForces()
        s.integrateEOM()
        s.update()
        print('equil:', i)
    return s


def calcDistances(coords, boxlength):
    distList = []
    nParticles = len(coords)
    for i in range(nParticles):
        for j in range(i+1, nParticles):
            distVec = [coords[i][0]-coords[j][0], coords[i][1]-coords[j][1]]
            distVec[0] -= boxlength*round(distVec[0]/boxlength)
            distVec[1] -= boxlength*round(distVec[1]/boxlength)
            dist = Norm(distVec)
            if dist < boxlength/2:
                distList.append(dist)
    return distList

def calcGr(coordList, boxlength, dr):
    distList = []
    density = (len(coordList[0][1])/(boxlength**2))
    for item in coordList:
        distList = distList+calcDistances(item[1], boxlength)
    nSample = 0.5*len(coordList[0][1])*len(coordList)
    gr=BinCounts(distList, dr)


    pairs = []
    for i in range(len(gr)):
        distList.append(Roundup(min(distList)-dr, dr)+(i-1)*dr)
    for i in range(len(gr)):
        gr[i] /= (nSample*density*math.pi*(((distList[i]+dr)**2)-(distList[i])**2))
        pairs.append([distList[i], gr[i]])
    return pairs





# Example of running the simulation and saving the data:

#PElist16, KElist16, TElist16, coordlist16 = runSimulation(100, 16, 2.5, 0.0005, 1.5, 500, 5000, 50, 100)
#print(calcGr(coordlist16, 16, 0.01))


# Example of creating an animated system that runs forever:
def genSystem(nparticles, boxlength, temperature, rcut, dt, colors, particleSize, bg, cop):
    s = System(nparticles, boxlength, temperature, rcut, dt, colors, particleSize, bg, cop)
    s.initVelocities()
    s.computeForces()
    s.integrateEOM()
    return s

def runAnimation(system):
    system.animate()
