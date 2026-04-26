Analytic study of the XY model The XY model is a system of two-component spins S2
x + S2
y = 1 with
an O(2) symmetry and nearest-neighbours interactions. The state of the spins can be characterized by a
single angle that describes the orientation of the vector (Sx, Sy). Study the behaviour of the system in two
dimensions by considering the high, and low temperature expansions. In the low temperature phase, the
fluctuations of the angles define the vorticity of the configuration. Study the contribution of the vortices to
the free energy of the system, and deduce an estimate of the critical temperature.
Introduce the Villain action for the XY model, and discuss the RG flow for this theory


−βF = log(Z(β))
S = − ∂F/∂T

P(X) = 1/Z * exp(-E/kT)
F = U - TS
S = kln(omega)

Run simulations by updating the angles in each grid point by a small dTheta x
Run the system in high temperature and low temperature configurations x
compute Vortex Density over time x
Graph the free energy of the system over time 
study the contribution of vortices to the free energy 
estimate the critical temperature 

write up results in a "short" report after doing some background reading.

    def applyFlipProbabilities(self):

        gridBefore = self.grid
        gridHamiltonianBefore = self.getAllSitesHamiltonian(gridBefore)
        
        gridAfter = self.grid + np.random.normal(loc=0.0, scale=15.0, size=(self.size,self.size))
        gridHamiltonianAfter =  self.getAllSitesHamiltonian(gridAfter)

        gridChangedHamiltonian = gridHamiltonianAfter - gridHamiltonianBefore
        gridChangeProbability  = np.exp((-1)*self.betaTemp*gridChangedHamiltonian)

        gridNew = np.zeros(shape=(self.size, self.size))

        for y in range(0, self.size):
            for x in range(0,self.size):

                if gridChangeProbability[y,x] >= np.random.random():
                    gridNew[y,x] = gridAfter[y,x] % 360
                else:
                    gridNew[y,x] = gridBefore[y,x] % 360

        self.grid = gridNew


dF = kT ln(Pa/Pb)
O(2) 

superfluid
helium, thin-films, superconductivity, liquid crystals, melting of two-dimensional crystals,
and dielectric plasma transition from a dielectric phase (charges bound as neutral dipoles)
to the conducting phase (free charges

    def getAnglePerturbation(self):
        return
np.random.normal(loc=0.0, scale=8.0)

    def applyFlipProbabilities(self):

        gridBefore = gridAfter = self.grid
        gridAfter += np.random.normal(loc=0.0, scale=8.0, size=(self.size, self.size))

        gridChangeProbability = self.getChangeProbablity(gridBefore, gridAfter)
        gridNew = np.zeros(shape=(self.size, self.size))

        for y in range(0, self.size):
            for x in range(0,self.size):
                if gridChangeProbability[y,x] >= np.random.random():
                    gridNew[y,x] = gridAfter[y,x] % 360

        self.grid = gridNew


    def getChangeEnergy(self, gridBefore, gridAfter):
        gridHamiltonianBefore  = self.getAllSitesHamiltonian(gridBefore)
        gridHamiltonianAfter   = self.getAllSitesHamiltonian(gridAfter)
        gridChangedHamiltonian = gridHamiltonianAfter - gridHamiltonianBefore
        return gridChangedHamiltonian

    def getChangeProbablity(self, gridBefore, gridAfter):
        gridChangedHamiltonian = self.getChangeEnergy(gridBefore, gridAfter)
        gridChangeProbability  = np.exp((-1)*self.betaTemp*gridChangedHamiltonian)
        return gridChangeProbability
       # plt.plot()
        # plt.figure(figsize=(16,9), dpi=160)
        # # plt.plot(self.betaTempTime, label="temp", linewidth=0.5)
        # plt.grid(linestyle=":", alpha=0.7)
        # plt.savefig(f"imgs/{fileName}_beta.png")

        # plt.close()
        # plt.figure(figsize=(16,9), dpi=160)
        # plt.xlabel("Beta")
        # plt.ylabel("Free Energy")
        # # plt.plot(self.betaTempTime, self.FreeEnergy, label="Free", linewidth=0.5)
        # plt.savefig(f"imgs/{fileName}_freeEnergy.png")
        # plt.close()

        # df = pd.DataFrame(columns={"BetaTemp":self.betaTempTime,"FreeEnergy":self.FreeEnergy})
        # df.to_csv(f"data/{fileName}.csv")