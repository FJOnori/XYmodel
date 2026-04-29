import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation # type: ignore

Kb = 1.380649*(10**(-23))

class Lattice():

    def __init__(self, size, betaTemp, couplingConstant, gridInit):
        self.size             = size
        self.betaTemp         = betaTemp
        self.couplingConstant = couplingConstant
        self.step             = 0

        if gridInit == "random":
            self.grid = np.array( [[float(np.random.randint(0,359)) \
                                   for _ in range(0, self.size)] \
                                   for _ in range(0,self.size)] )
        elif gridInit == "uniform":
            angle = float(np.random.randint(0,359))
            self.grid = np.array( [[ angle for _ in range(0, self.size)] \
                                   for _ in range(0,self.size)] )

        self.VortexNumberFrames    = []
        self.VortexAnimationFrames = []
        self.NumberFrames          = []
        self.Animationframes       = []
        self.FreeEnergy            = []
        self.betaTempTime          = []
        self.vortexDensity         = []
        self.vortexNetDensity      = []
    
    def getSiteKernel(self, coordinate, inputGrid=None):
        
        d  = self.size
        x  = int((coordinate[0]) % d)
        y  = int((coordinate[1]) % d)
        xL = int((coordinate[0] - 1) % d)
        xR = int((coordinate[0] + 1) % d)
        yU = int((coordinate[1] - 1) % d)
        yD = int((coordinate[1] + 1) % d)

        if inputGrid is None: lattice = self.grid
        else: lattice = inputGrid
        
        kernel = np.array([[ lattice[yU,xL], lattice[yU,x], lattice[yU,xR] ],
                           [ lattice[y, xL], lattice[y, x], lattice[y, xR] ],
                           [ lattice[yD,xL], lattice[yD,x], lattice[yD,xR] ]])
        
        return kernel
        
    def getSiteHamiltonian(self, coordinate, inputGrid=None):

        if inputGrid is None: lattice = self.grid
        else: lattice = inputGrid
        
        d  = self.size
        x  = int((coordinate[0]) % d)
        y  = int((coordinate[1]) % d)

        kernel = self.getSiteKernel(coordinate, inputGrid)

        kernelDiff = np.cos( (kernel - lattice[y, x]) * (np.pi/180) )
        siteHamiltonian = (-self.couplingConstant)*(np.sum(kernelDiff) - 1) 

        return siteHamiltonian
    
    def getAllSitesHamiltonian(self, inputGrid=None):

        if inputGrid is None: lattice = self.grid
        else: lattice = inputGrid

        gridHamiltonian = np.array( [[0.0 for _ in range(0, self.size)] for _ in range(0,self.size)] )

        for y in range(0, self.size):
            for x in range(0,self.size):
                gridHamiltonian[y,x] += self.getSiteHamiltonian(coordinate=(x,y), inputGrid=lattice)

        return gridHamiltonian

    def getFullHamiltonian(self,inputGrid=None):
        return (np.sum(self.getAllSitesHamiltonian(inputGrid))/2)
    
    def getAnglePerturbation(self):
        return np.random.normal(loc=0.0, scale=8.0)

    def getSingleFlipProbability(self, coordinate, perturbation=None, inputGrid=None):

        if inputGrid is None: lattice = self.grid
        else: lattice = inputGrid

        if perturbation is None: dTheta = np.random.normal(loc=0.0, scale=8.0)
        else: dTheta = perturbation

        x  = int((coordinate[0]) % self.size)
        y  = int((coordinate[1]) % self.size)

        gridBefore = lattice.copy()
        gridAfter  = lattice.copy()

        gridAfter[y,x] += dTheta

        HamiltonianBefore  = self.getSiteHamiltonian(coordinate, inputGrid=gridBefore)
        HamiltonianAfter   = self.getSiteHamiltonian(coordinate, inputGrid=gridAfter)

        ChangedHamiltonian = HamiltonianBefore - HamiltonianAfter
        q = np.exp(self.betaTemp*ChangedHamiltonian)
        ChangeProbability  = q / (1+q)

        return ChangeProbability, gridAfter

    def applySingleFlipProbability(self, coordinate, perturbation=None, inputGrid=None):
        
        if inputGrid is None: lattice = self.grid
        else: lattice = inputGrid

        if perturbation is None: dTheta = np.random.normal(loc=0.0, scale=8.0)
        else: dTheta = perturbation

        ChangeProbability, gridAfter = self.getSingleFlipProbability(coordinate, perturbation=dTheta, inputGrid=lattice)
        
        x  = int((coordinate[0]) % self.size)
        y  = int((coordinate[1]) % self.size)

        if ChangeProbability >= np.random.random():
            self.grid[y,x] = gridAfter[y,x] % 360

    def getUnitVector(self, angle):
        angleDegrees = angle * (np.pi/180)
        return np.array([np.cos(angleDegrees), np.sin(angleDegrees)])

    def angleDiff(self, theta2, theta1):
        d = theta2 - theta1
        return (d + np.pi) % (2*np.pi) - np.pi

    def getSiteVortex(self, coordinate, inputGrid = None):

        if inputGrid is None: lattice = self.grid
        else: lattice = inputGrid
    
        d = self.size
        lattice = self.grid * (np.pi / 180)  # convert to radians

        # periodic boundaries
        x  = int((coordinate[0]) % d)
        y  = int((coordinate[1]) % d)
        xR = int((coordinate[0] + 1) % d)
        yD = int((coordinate[1] + 1) % d)

        t1 = lattice[y,  x]
        t2 = lattice[y,  xR]
        t3 = lattice[yD, xR]
        t4 = lattice[yD, x]

        winding = (
            self.angleDiff(t2, t1) +
            self.angleDiff(t3, t2) +
            self.angleDiff(t4, t3) +
            self.angleDiff(t1, t4)
        )

        return winding / (2 * np.pi)

    def getAllSitesVortex(self, inputGrid=None):

        if inputGrid is None: lattice = self.grid
        else: lattice = inputGrid

        gridVortex = np.array( [[0.0 for _ in range(0, self.size)] for _ in range(0,self.size)] )

        for y in range(0, self.size):
            for x in range(0,self.size):

                v = self.getSiteVortex(coordinate=[x,y], inputGrid=lattice)
                
                if v > 0.5:
                    gridVortex[y,x] = +1   # vortex
                elif v < -0.5:
                    gridVortex[y,x] = -1   # antivortex
                else:
                    gridVortex[y,x] =  0
                
        return gridVortex
        
    def getFullVortexDensity(self, inputGrid=None):
        count = 0

        for y in range(self.size):
            for x in range(self.size):
                coordinate = [x,y]
                v = self.getSiteVortex(coordinate)
                if abs(v) > 0.5:
                    count += 1

        return count / (self.size**2)

    def runSimulation(self, steps):
        
        self.step = 0
        self.NumberFrames = []
        self.VortexNumberFrames = []

        for n in range(0,steps):
            self.step += 1
            print(f"Running Simulation: {np.round(100*n/steps,3):.3f} %", end="\r")

            for i in range(0, self.size**2):
                coordinate = [int(np.random.randint(0,self.size)), int(np.random.randint(0,self.size))]
                self.applySingleFlipProbability(coordinate=coordinate)
            

            self.step += 1
            # self.betaTemp += 0.0005
            self.betaTempTime.append(self.betaTemp)
            self.vortexDensity.append( np.sum(np.abs(self.getAllSitesVortex())) / (2 * (self.size**2)) )
            self.vortexNetDensity.append(np.sum(self.getAllSitesVortex()))

            print(f"Running Simulation: 100.00 %", end="\r")
            self.NumberFrames.append(self.grid.copy())
            self.VortexNumberFrames.append(self.getAllSitesVortex().copy())
        
    def createImage(self, frame):
        plt.xticks([])
        plt.yticks([])
        im = plt.imshow(self.NumberFrames[frame], cmap="hsv", animated=True)
        txt = plt.text(-0.5, -0.9, f"β = {np.round(self.betaTempTime[frame],4):.4f}\nstep = {frame}")
        return [im,txt]

    def createVortexImage(self, frame):
        plt.xticks([])
        plt.yticks([])
        im = plt.imshow(self.VortexNumberFrames[frame], cmap="bwr", animated=True)
        txt = plt.text(-0.5, -0.9, f"β = {np.round(self.betaTempTime[frame],4):.4f}\nstep = {frame}")
        return [im,txt]

    def createAnimation(self, fileName, aniType):
        
        fig, ax = plt.subplots()
        ax.set_xticks([])
        ax.set_yticks([])

        self.Animationframes = []
        steps = len(self.NumberFrames)

        print()
        for frame in range(0,steps):
            print(f"{fileName} Rendering Frames: {np.round(100*frame/steps,2):.2f} %", end="\r")
            if aniType == "phase":
                self.Animationframes.append(self.createImage(frame).copy())
            elif aniType == "vortex":
                self.Animationframes.append(self.createVortexImage(frame).copy())

        print(f"{fileName} Rendering Frames: 100.00 %", end="\r")
        print()
        print("Creating animation...")
        ani = animation.ArtistAnimation(fig, self.Animationframes, interval=60, blit=True, repeat_delay=1000)
        ani.save(f"videos/{fileName}.mp4", writer = "ffmpeg", fps = 24)
        print("Done.")

    def saveObservables(self, fileName):
        a = np.array([self.betaTempTime, self.vortexDensity, self.vortexNetDensity])
        np.savetxt(f"data/{fileName}.csv", a, delimiter=",")


if __name__ == "__main__":
    sim = Lattice(size=50, betaTemp=1/0.893, couplingConstant=1, gridInit="random")
    sim.runSimulation(steps=1000)
    sim.saveObservables("crit50")
    sim.createAnimation(fileName="crit50", aniType="phase")
    sim.createAnimation(fileName="crit50Vortex", aniType="vortex")
