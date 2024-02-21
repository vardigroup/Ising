from os import write
import numpy as np
import random
from boolean_formula import Formula
from itertools import product

class IsingModel:
    def __init__(self, beta = 1, mu = 1, interactions = []):
        self._numLatticeSites = len(interactions)
        self._beta = beta              # Inverse temperature
        self._mu = mu                # Field orientation
        self._interactions = interactions     # A square upper triangular matrix - diagonal entries are h, others are J

    def toWMC(self):
        """
        Creates a weighted model counting problem whose solution
        is the partition function of the Ising model.

        :return: A Formula which represents the weighted model counting
        problem 
        """
        form = Formula()
        varIds = [[0 for i in range(self._numLatticeSites)] for i in range(self._numLatticeSites)]
        
        # Create variables for each lattice site and each non-zero pairwise interaction between them
        for i in range(self._numLatticeSites):
            varIds[i][i] = form.fresh_variable(np.exp(-1 * self._mu * self._beta * self._interactions[i][i]), np.exp(self._mu * self._beta * self._interactions[i][i]))
            for j in range(i + 1, self._numLatticeSites):
                if self._interactions[i][j] != 0:
                    varIds[i][j] = form.fresh_variable(np.exp(-1 * self._beta * self._interactions[i][j]), np.exp(self._beta * self._interactions[i][j]))
        
        # Add clauses for each pairwise interaction
        for i in range(self._numLatticeSites):
            for j in range(i + 1, self._numLatticeSites):
                if varIds[i][j] != 0:
                    form.add_clause([varIds[i][j], varIds[i][i], varIds[j][j]])
                    form.add_clause([varIds[i][j], -varIds[i][i], -varIds[j][j]])
                    form.add_clause([-varIds[i][j], varIds[i][i], -varIds[j][j]])
                    form.add_clause([-varIds[i][j], -varIds[i][i], varIds[j][j]])
        return form

    def to_UAI08(self, filename):
        out = open(filename, "w")

        ### PREAMBLE
        # First part of preamble
        out.write("ISING\n")
        out.write(str(self._numLatticeSites) + "\n")
        out.write("2 " * (self._numLatticeSites - 1) + "2\n")
        out.write(str(self.numUnaryFuncs()) + " " + str(self.numBinaryFuncs()) + " " + str(self._beta) + " " + str(self._mu) + "\n")

        # Write function inputs part of preamble 
        for i in range(self._numLatticeSites):
            if self._interactions[i][i] != 0:
                out.write("1 " + str(i) + "\n")
        for i in range(self._numLatticeSites):
            for j in range(i + 1, self._numLatticeSites):
                if self._interactions[i][j] != 0:
                    out.write("2 " + str(i) + " " + str(j) + "\n")

        out.write("\n")
        # FUNCTION TABLE
        for i in range(self._numLatticeSites):
            if self._interactions[i][i] != 0:
                out.write("2\n")
                out.write(" " + str(-self._interactions[i][i]) + " " + str(self._interactions[i][i]) + "\n")
        for i in range(self._numLatticeSites):
            for j in range(i + 1, self._numLatticeSites):
                if self._interactions[i][j] != 0:
                    out.write("4\n")
                    out.write(" " + str(self._interactions[i][j]) + " " + str(-self._interactions[i][j]) + "\n")
                    out.write(" " + str(-self._interactions[i][j]) + " " + str(self._interactions[i][j]) + "\n")
        out.close()

    """
    Outputs the Ising model as a graph in the two-file representation used by
    https://github.com/panzhang83/catn
    """
    def to_pan_format(self, fileprefix):
        # Writes nodes file to describe edges present in model
        with open('{}nodes.txt'.format(fileprefix), "w") as graph:
            graph.write(str(self._numLatticeSites) + " " + str(self.numBinaryFuncs()) + "\n")
            for i in range(self._numLatticeSites):
                for j in range(i + 1, self._numLatticeSites):
                    if self._interactions[i][j] != 0:
                        graph.write(str(i) + " " + str(j) + "\n")
                            
        # Writes Jij file to give weights of each interaction
        Jvals = np.matrix(np.array(self._interactions))
        with open('Jij{}nodes.txt'.format(fileprefix), "w") as Jnodes:
            for line in Jvals:
                np.savetxt(Jnodes, line, fmt='%.2f')
         
                  
    '''
    Same behavior as readline, but discards lines starting with a #
    and parts of lines after a #
    '''
    @staticmethod
    def readline_comment(input):
        line = ""
        while line == "":
            line = input.readline().split("#")[0]
        return line


    @staticmethod
    def from_UAI08(in_model):
        ### PREAMBLE
        # Read meta-information
        if not IsingModel.readline_comment(in_model).startswith("ISING"):
            raise Exception("Misformated file " + in_model)
        numLatticeSites = int(IsingModel.readline_comment(in_model))
        IsingModel.readline_comment(in_model)         # Disregard line specifying variable domains
        # The next line may contain optional fields
        metaline = [int(float(i)) for i in IsingModel.readline_comment(in_model).split()]
        num_h_vals, num_J_vals = metaline[:2]
        beta = 1
        mu = 1
        if len(metaline) > 2:
            beta = float(metaline[2])
            if len(metaline) > 3:
                mu = int(metaline[3])
                
        # Function inputs
        nonzero_h_slots = []
        nonzero_J_slots = []
        for i in range(num_h_vals):
            nonzero_h_slots.append(int(IsingModel.readline_comment(in_model).split()[1]))
        for i in range(num_J_vals):
            nonzero_J_slots.append(tuple([int(j) for j in IsingModel.readline_comment(in_model).split()[1:]]))


        IsingModel.readline_comment(in_model)         # Read empty line
        ### FUNCTION TABLE
        interactions = [[0 for i in range(numLatticeSites)] for j in range(numLatticeSites)]
        # H values
        for i in nonzero_h_slots:
            if IsingModel.readline_comment(in_model) != "2\n":
                #error
                pass
            interactions[i][i] = float(IsingModel.readline_comment(in_model).split()[1])
        # J values
        for (i,j) in nonzero_J_slots:
            if IsingModel.readline_comment(in_model) != "4\n":
                #error
                pass
            interactions[i][j] = float(IsingModel.readline_comment(in_model).split()[0])
            IsingModel.readline_comment(in_model)
        
        in_model.close()
        return(IsingModel(beta = beta, mu = mu, interactions = interactions))

    def numUnaryFuncs(self):
        unaryCount = 0
        for i in range(self._numLatticeSites):
            if self._interactions[i][i] != 0:
                unaryCount += 1
        return unaryCount

    def numBinaryFuncs(self):
        numBinaryFuncs = 0
        for i in range(self._numLatticeSites):
            for j in range(i + 1, self._numLatticeSites):
                if self._interactions[i][j] != 0:
                    numBinaryFuncs += 1
        return numBinaryFuncs
        
    """
    Generates an Ising model comprised of an x by y by z grid of lattice sites with 
    edge J values given by J(x1,y1,z1,x2,y2,z2)->float and h=0

    f should be a symetric function giving the full strength of the interaction between the inputs
    """        
    @staticmethod
    def ThreeDGrid(x, y, z, f, beta):
        model = IsingModel(beta = beta, interactions= [[0 for j in range(x*y*z)] for i in range(x*y*z)])
        for (i1,j1,k1) in product(range(x),range(y),range(z)):
            for(i2,j2,k2) in product(range(x),range(y),range(z)):
                if (i1, j1, k1) == (i2, j2, k2):
                    continue
                model._interactions[i1*y*z + j1*z + k1][i2*y*z + j2*z + k2] = f(i1,j1,k1,i2,j2,k2)
        
        # Correct for double counting
        for i in range(len(model._interactions)):
            for j in range(i+1, len(model._interactions)):
                model._interactions[j][i] = 0

        return model

    """
    Generates an Ising model comprised of an x by y grid of lattice sites with 
    edge J values given by J(x1,y1,x2,y2)->float and h=0

    f should be a symmetric function giving the full strength of the interaction between the inputs
    """        
    @staticmethod
    def TwoDGrid(x, y, f, beta):
        model = IsingModel(beta = beta, interactions= [[0 for j in range(x*y)] for i in range(x*y)])
        for (i1,j1) in product(range(x),range(y)):
            for(i2,j2) in product(range(x),range(y)):
                if (i1, j1) == (i2, j2):
                    continue
                model._interactions[i1*y + j1][i2*y + j2] = f(i1,j1,i2,j2)
        
        # Correct for double counting
        for i in range(len(model._interactions)):
            for j in range(i+1, len(model._interactions)):
                model._interactions[j][i] = 0

        return model

    """
    Generates an Ising model as a random graph of n nodes with 
    expected degree expDegree using the seed value seed.
    The upper and lower bounds on J are given by JLB and JUB respectively, and J values
    are sampled uniformly between them.

    Does not include field (h = 0).
    """        
    @staticmethod
    def random(n: int, expDegree: float, beta: float, JLB: float = -1, JUB: float = 1, seed: int = 0):
        random.seed(seed)
        model = IsingModel(beta = beta, interactions= [[0 for j in range(n)] for i in range(n)])
        for j in range(n):
            for i in range(j):
                model._interactions[i][j] = random.choices([0,1], weights=(n - expDegree, expDegree), k=1)[0] * random.uniform(JLB, JUB)
        print(model._interactions)

        return model

    '''
    Produces a formula describing a 2D Ising model where, instead of having constraints relating vertices to edges, we have clauses for each face/square in the grid.
    One can check that this encoding is equivalent. Hopefully, reducing the number of variables will make stuff run faster. 
    '''
    @staticmethod
    def TwoDCondense(n : int):
        form = Formula()
        varIdsHoriz = [[form.fresh_variable(np.exp(-1), np.exp(1)) for j in range(n)] for i in range(n-1)]
        varIdsVert = [[form.fresh_variable(np.exp(-1), np.exp(1)) for j in range(n-1)] for i in range(n)]
        
        # Add clauses for each face of the grid
        for i in range(n-1):
            for j in range(n-1):
                form.add_clause([-varIdsHoriz[i][j], varIdsHoriz[i][j+1], varIdsVert[i][j], varIdsVert[i+1][j]])
                form.add_clause([varIdsHoriz[i][j], -varIdsHoriz[i][j+1], varIdsVert[i][j], varIdsVert[i+1][j]])
                form.add_clause([varIdsHoriz[i][j], varIdsHoriz[i][j+1], -varIdsVert[i][j], varIdsVert[i+1][j]])
                form.add_clause([varIdsHoriz[i][j], varIdsHoriz[i][j+1], varIdsVert[i][j], -varIdsVert[i+1][j]])
                form.add_clause([varIdsHoriz[i][j], -varIdsHoriz[i][j+1], -varIdsVert[i][j], -varIdsVert[i+1][j]])
                form.add_clause([-varIdsHoriz[i][j], varIdsHoriz[i][j+1], -varIdsVert[i][j], -varIdsVert[i+1][j]])
                form.add_clause([-varIdsHoriz[i][j], -varIdsHoriz[i][j+1], varIdsVert[i][j], -varIdsVert[i+1][j]])
                form.add_clause([-varIdsHoriz[i][j], -varIdsHoriz[i][j+1], -varIdsVert[i][j], varIdsVert[i+1][j]])
        return form