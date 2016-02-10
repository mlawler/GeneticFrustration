from pylab import *
import time

N1 = 10
N2 = 10
Npop = 100
mutation_rate=0.01 #careful with this: a high mutation rate can saturate the population with mutants!
Ngen=100
SpeciesCorr=0.1 # A correlation coefficient greater than this is significant
SpeciesPop=20
lamfrac=0.1


########## Definition of species ##################
# a species is an object! It has well defined rules for what it can
# do and what its properties are.

# Note: __X is how you declare a private variable. Can we declare private functions too?

# J1-J2-J3-J3p Heisenberg model
symJ=array([1.0,0.2,0.01,0.01]) # The DNA for a particular J1-J2-J3-J3p model
asymJ=ones(18+18+9+18)
asymJ[18:36]=0.2*ones(18)
asymJ[36:45]=0.01*ones(9)
asymJ[45:]=0.01*ones(18)

########## Structures of the Hamiltonian ##############
# These should be moved to a file (eventually)

symTerms=[
    [ # J1 term
        [0,3,0,0],[1,4,0,0],[2,5,0,0],[0,6,0,0],[1,7,0,0],[2,8,0,0],
        [3,0,0,0],[4,1,0,0],[5,2,0,0],[3,6,0,0],[4,7,0,0],[5,8,0,0],
        [6,0,0,0],[7,1,0,0],[8,2,0,0],[6,3,0,0],[7,4,0,0],[8,5,0,0],
        [3,0,1,0],[4,1,1,0],[5,2,1,0],[0,3,-1,0],[1,4,-1,0],[2,5,-1,0],
        [6,3,0,1],[7,4,0,1],[8,5,0,1],[3,6,0,-1],[4,7,0,-1],[5,8,0,-1],
        [6,0,1,1],[7,1,1,1],[8,2,1,1],[0,6,-1,-1],[1,7,-1,-1],[2,8,-1,-1]
    ],[ # J2 term
        [6,0,1,0],[7,1,1,0],[8,2,1,0],[3,6,1,0],[4,7,1,0],[5,8,1,0],
        [0,6,-1,0],[1,7,-1,0],[2,8,-1,0],[6,3,-1,0],[7,4,-1,0],[8,5,-1,0],        
        [0,3,0,1],[1,4,0,1],[2,5,0,1],[6,0,0,1],[7,1,0,1],[8,2,0,1],
        [0,6,0,-1],[1,7,0,-1],[2,8,0,-1],[3,0,0,-1],[4,1,0,-1],[5,2,0,-1],
        [3,0,1,1],[4,1,1,1],[5,2,1,1],[6,3,1,1],[7,4,1,1],[8,5,1,1],
        [0,3,-1,-1],[1,4,-1,-1],[2,5,-1,-1],[3,6,-1,-1],[4,7,-1,-1],[5,8,-1,-1]
    ],[ # J3 terms
        [6,6,1,0],[7,7,1,0],[8,8,1,0],[6,6,-1,0],[7,7,-1,0],[8,8,-1,0],
        [0,0,0,1],[1,1,0,1],[2,2,0,1],[0,0,0,-1],[1,1,0,-1],[2,2,0,-1],
        [3,3,1,1],[4,4,1,1],[5,5,1,1],[3,3,-1,-1],[4,4,-1,-1],[5,5,-1,-1]
    ],[ # J3p terms
        [0,0,1,0],[1,1,1,0],[2,2,1,0],[3,3,1,0],[4,4,1,0],[5,5,1,0],
        [0,0,-1,0],[1,1,-1,0],[2,2,-1,0],[3,3,-1,0],[4,4,-1,0],[5,5,-1,0],
        [3,3,0,1],[4,4,0,1],[5,5,0,1],[6,6,0,1],[7,7,0,1],[8,8,0,1],
        [3,3,0,-1],[4,4,0,-1],[5,5,0,-1],[6,6,0,-1],[7,7,0,-1],[8,8,0,-1],
        [0,0,1,1],[1,1,1,1],[2,2,1,1],[6,6,1,1],[7,7,1,1],[8,8,1,1],
        [0,0,-1,-1],[1,1,-1,-1],[2,2,-1,-1],[6,6,-1,-1],[7,7,-1,-1],[8,8,-1,-1],
    ]
]

asymTerms=[ # J1 terms
    [[0,3,0,0],[3,0,0,0]],
    [[1,4,0,0],[4,1,0,0]],
    [[2,5,0,0],[5,2,0,0]],
    [[0,6,0,0],[6,0,0,0]],
    [[1,7,0,0],[7,1,0,0]],
    [[2,8,0,0],[8,2,0,0]],
    [[3,6,0,0],[6,3,0,0]],
    [[4,7,0,0],[7,4,0,0]],
    [[5,8,0,0],[8,5,0,0]],
    [[3,0,1,0],[0,3,-1,0]],
    [[4,1,1,0],[1,4,-1,0]],
    [[5,2,1,0],[2,5,-1,0]],
    [[6,3,0,1],[3,6,0,-1]],
    [[7,4,0,1],[4,7,0,-1]],
    [[8,5,0,1],[5,8,0,-1]],
    [[6,0,1,1],[0,6,-1,-1]],
    [[7,1,1,1],[1,7,-1,-1]],
    [[8,2,1,1],[2,8,-1,-1]],
    # J2 terms
    [[6,0,1,0],[0,6,-1,0]],
    [[7,1,1,0],[1,7,-1,0]],
    [[8,2,1,0],[2,8,-1,0]],
    [[3,6,1,0],[6,3,-1,0]],
    [[4,7,1,0],[7,4,-1,0]],
    [[5,8,1,0],[8,5,-1,0]],        
    [[0,3,0,1],[3,0,0,-1]],
    [[1,4,0,1],[4,1,0,-1]],
    [[2,5,0,1],[5,2,0,-1]],
    [[6,0,0,1],[0,6,0,-1]],
    [[7,1,0,1],[1,7,0,-1]],
    [[8,2,0,1],[2,8,0,-1]],
    [[3,0,1,1],[0,3,-1,-1]],
    [[4,1,1,1],[1,4,-1,-1]],
    [[5,2,1,1],[2,5,-1,-1]],
    [[6,3,1,1],[3,6,-1,-1]],
    [[7,4,1,1],[4,7,-1,-1]],
    [[8,5,1,1],[5,8,-1,-1]],
    # J3 terms
    [[6,6,1,0],[6,6,-1,0]],
    [[7,7,1,0],[7,7,-1,0]],
    [[8,8,1,0],[8,8,-1,0]],
    [[0,0,0,1],[0,0,0,-1]],
    [[1,1,0,1],[1,1,0,-1]],
    [[2,2,0,1],[2,2,0,-1]],
    [[3,3,1,1],[3,3,-1,-1]],
    [[4,4,1,1],[4,4,-1,-1]],
    [[5,5,1,1],[5,5,-1,-1]],
    # J3p terms
    [[0,0,1,0],[0,0,-1,0]],
    [[1,1,1,0],[1,1,-1,0]],
    [[2,2,1,0],[2,2,-1,0]],
    [[3,3,1,0],[3,3,-1,0]],
    [[4,4,1,0],[4,4,-1,0]],
    [[5,5,1,0],[5,5,-1,0]],
    [[3,3,0,1],[3,3,0,-1]],
    [[4,4,0,1],[4,4,0,-1]],
    [[5,5,0,1],[5,5,0,-1]],
    [[6,6,0,1],[6,6,0,-1]],
    [[7,7,0,1],[7,7,0,-1]],
    [[8,8,0,1],[8,8,0,-1]],
    [[0,0,1,1],[0,0,-1,-1]],
    [[1,1,1,1],[1,1,-1,-1]],
    [[2,2,1,1],[2,2,-1,-1]],
    [[6,6,1,1],[6,6,-1,-1]],
    [[7,7,1,1],[7,7,-1,-1]],
    [[8,8,1,1],[8,8,-1,-1]],
]

class Ham:
    def __init__(self, oldHam=None, model=[], genome=array([])):
        if oldHam:
            self.model=copy(oldHam.model)
            self.genome=oldHam.genome.copy()
        elif model:
            mlen=len(model)
            self.model = model
            self.genome=rand(mlen)-0.5*ones(mlen)
        else:
            self.model = model
            self.genome = genome

    def FourierTransform(self):
        tildeJ = zeros((9,9,N1,N2),dtype=complex)
    
        for i in range(len(self.model)):
            for term in self.model[i]:
                a,b,m1,m2 = term
                for n1 in range(N1):
                    for n2 in range(N2):
                        tildeJ[a,b,n1,n2] = self.genome[i]*exp(-1.0j*2*pi*(m1*n1+m2*n2)/(N1*N2))
        return tildeJ

    def LowDOS(self):
        tildeJ = FourierTransform(self.J)
    
        lams = []
        for n1 in range(N1):
            for n2 in range(N2):
                lams.extend(eigvalsh(tildeJ[:,:,n1,n2]))
    
        lams.sort()
        lammin = lams[0]
        lammax = lams[-1]
        W = lammax-lammin  # Bandwidth
    
        dos = 0
        for lam in lams:
            if lam-lammin < lamfrac*W: #lamfrac is a global variable here
                dos = dos + 1

        return dos #/ (1.0*9*N1*N2)        

    
############ Initialization ####################

def Initialization(microH=None):
    igen=[]
    for i in range(Npop):
        if microH:
            igen.append(Ham(microHam=microH))
        else:
            igen.append(Ham())
    return igen
   
############# Fitness ####################

def Fitness(gen,microHam):
    return IndividualFitness(gen,microHam)

def IndividualFitness(gen,microHam):
    fitlist = zeros(Npop)
    for i in range(Npop):
        dist=Distance(gen[i],microHam)
        lDOS = gen[i].LowDOS()
        fitlist[i] = lDOS/(1.0+dist) # Ensure that dist < 1.0 is close enough for search
    return fitlist

def LowResourcesFitness(gen,microHam):
    # This seems to be very sloooooow
    fitlist = zeros(Npop)
    for i in range(Npop):
        Nspecies=0
        for j in range(i+1,Npop):
            if Correlation(gen[i],gen[j]) > SpeciesCorr:
                Nspecies= Nspecies + 1
        dist = Distance(gen[i],microHam)
        lDOS = gen[i].LowDOS()
        resources = exp(-Nspecies**2/(2.0*SpeciesPop**2))
        fitlist[i] = lDOS/(1.0+dist)*resources
    return fitlist

def Correlation(H1,H2):
    fJ1=list(H1.J.flatten())
    fJ2=list(H2.J.flatten())
    fK1=[]
    fK2=[]
    for i in range(len(fJ1)):
        if fJ1[i] != 0 or fJ2[i] != 0:
            fK1.append(fJ1[i])
            fK2.append(fJ2[i])
    return corrcoef(fK1,fK2)[0,1]

def Distance(H1,H2):
    AnticorrDist = 2.0 # The distance of a perfectly negative correlation
    corr = Correlation(H1,H2)
    if corr > 0.0:
        return 1/corr-1
    elif corr < 0.0:
        return -AnticorrDist/corr
    else:
        return inf

def GroupSpecies(gen):
    groups = []
    ungrouped = list(gen)
    
    grouped=[]
    
    def FindGroup(ungrouped):
        group=[]
        remaining=[]
        example=ungrouped.pop()
        group.append(example)
        while len(ungrouped)>0:
            test=ungrouped.pop()
            if Correlation(example,test) > SpeciesCorr:
                group.append(test)
            else:
                remaining.append(test)
        return group, remaining
    
    while len(ungrouped) > 0:
        group,remaining = FindGroup(ungrouped)
        grouped.append(group)
        ungrouped=remaining
        
    return grouped

    
def Psig(N):
    # define two arrays as correlated if there is a 1 standard deviation chance.
    c=zeros(10000)
    for i in range(10000):
        c[i]=Correlation(rand(N),rand(N))
    return mean(c), std(c)

################ Selection ###################

def Selection(fitlist):
    return RouletteSelection(fitlist)

def RouletteSelection(fitlist):
    flcumsum=fitlist.cumsum()
    p=rand()*flcumsum[-1]
    for i,total in enumerate(flcumsum):
        if p<total:
            return i

def ThresholdSelection(fitlist):
    fave = average(fitlist)
    fstd = std(fitlist)

    # keep randomly selecting an individual until they are fit!
    i=random_integers(len(fitlist))-1
    while (fitlist[i]<fave+1.2*fstd):
        i = random_integers(len(fitlist))-1

    return i

################### Crossover ###################
# This preserves translational symmetry by symmetrizing the domains
# that are crossed over.

def Crossover(H1,H2):
    # H1 and H2 are the parents. The children are
    K1=Ham(H1)
    K2=Ham(H2)
    
    # Divide the rank 4 tensors H1.J and H2.J into 4*9=36 regions defined by randomly chosen planes:
    a=random_integers(8)
    m1=random_integers(1,N1/2)
    m2=random_integers(1,N2/2)

    # and then swap every other region
    K1.J[:a,a:,:m1,:m2]=H2.J[:a,a:,:m1,:m2]
    K1.J[:a,a:,(-m1%N1+1):,:m2]=H2.J[:a,a:,(-m1%N1+1):,:m2]
    K1.J[:a,a:,:m1:,(-m2%N2+1):]=H2.J[:a,a:,:m1,(-m2%N2+1):]
    K1.J[:a,a:,(-m1%N1+1):,(-m2%N2+1):]=H2.J[:a,a:,(-m1%N1+1):,(-m2%N2+1):]
    K1.J[:a,a:,m1:(-m1%N1+1),m2:(-m2%N2+1)]=H2.J[:a,a:,m1:(-m1%N1+1),m2:(-m2%N2+1)]

    K1.J[a:,:a,:m1,:m2]=H2.J[a:,:a,:m1,:m2]
    K1.J[a:,:a,(-m1%N1+1):,:m2]=H2.J[a:,:a,(-m1%N1+1):,:m2]
    K1.J[a:,:a,:m1:,(-m2%N2+1):]=H2.J[a:,:a,:m1,(-m2%N2+1):]
    K1.J[a:,:a,(-m1%N1+1):,(-m2%N2+1):]=H2.J[a:,:a,(-m1%N1+1):,(-m2%N2+1):]
    K1.J[a:,:a,m1:(-m1%N1+1),m2:(-m2%N2+1)]=H2.J[a:,:a,m1:(-m1%N1+1),m2:(-m2%N2+1)]
    
    K2.J[:a,a:,:m1,:m2]=H1.J[:a,a:,:m1,:m2]
    K2.J[:a,a:,(-m1%N1+1):,:m2]=H1.J[:a,a:,(-m1%N1+1):,:m2]
    K2.J[:a,a:,:m1:,(-m2%N2+1):]=H1.J[:a,a:,:m1,(-m2%N2+1):]
    K2.J[:a,a:,(-m1%N1+1):,(-m2%N2+1):]=H1.J[:a,a:,(-m1%N1+1):,(-m2%N2+1):]
    K2.J[:a,a:,m1:(-m1%N1+1),m2:(-m2%N2+1)]=H1.J[:a,a:,m1:(-m1%N1+1),m2:(-m2%N2+1)]

    K2.J[a:,:a,:m1,:m2]=H1.J[a:,:a,:m1,:m2]
    K2.J[a:,:a,(-m1%N1+1):,:m2]=H1.J[a:,:a,(-m1%N1+1):,:m2]
    K2.J[a:,:a,:m1:,(-m2%N2+1):]=H1.J[a:,:a,:m1,(-m2%N2+1):]
    K2.J[a:,:a,(-m1%N1+1):,(-m2%N2+1):]=H1.J[a:,:a,(-m1%N1+1):,(-m2%N2+1):]
    K2.J[a:,:a,m1:(-m1%N1+1),m2:(-m2%N2+1)]=H1.J[a:,:a,m1:(-m1%N1+1),m2:(-m2%N2+1)]
    
    return K1,K2

################## Mutation ####################
# This violates symmetry!!!!

def Mutation(h):
    return ProportionalFlipMutation(h)
    
def FlipMutation(h):
# randomly add a random number to any Jij
    K=Ham(h)
    for a in range(9):
        for b in range(9):
            for n1 in range(N1):
                for n2 in range(N2):
                    if rand() < mutation_rate:
                        K.J[a,b,n1,n2]=rand()-0.5
                        K.J[b,a,-n1%N1,-n2%N2]=K.J[a,b,n1,n2]
    return K

def ProportionalFlipMutation(h):
# randomly add a proportional random number to any Jij
    K=Ham(h)
    for a in range(9):
        for b in range(9):
            for n1 in range(N1):
                for n2 in range(N2):
                    if rand() < mutation_rate:
                        K.J[a,b,n1,n2]=h.J[a,b,n1,n2]*(rand()-0.5)
                        K.J[b,a,-n1%N1,-n2%N2]=K.J[a,b,n1,n2]
    return K

def SmallMutation(h):
    Jmax = h.J.max()
    K=Ham(h)
    for a in range(9):
        for b in range(9):
            for n1 in range(N1):
                for n2 in range(N2):
                    if rand() < mutation_rate:
                        K.J[a,b,n1,n2]=h.J[a,b,n1,n2]+(Jmax/10.0)*(rand()-0.5)
                        K.J[b,a,-n1%N1,-n2%N2]=K.J[a,b,n1,n2]
    return K
 


############## Main genetic alogirthm ##################   

def GeneticAlgorithm(gen,microHam,stats=[]):
    # Recording the fitness of the previous generation seems confusing.
    for g in range(Ngen):
        fitlist=Fitness(gen,microHam)

        nextgen = []
        for i in range(Npop/2):
            parent1 = Selection(fitlist)
            parent2 = Selection(fitlist)
            K1,K2 = Crossover(gen[parent1],gen[parent2])
#            print i, parent1, parent2, "Test H: ", TestHermetian(K1), TestHermetian(K2)
            child1 = Mutation(K1)
            child2 = Mutation(K2)
#            print "Test Hermetian: ", TestHermetian(child1), TestHermetian(child2)
            nextgen.append(child1)
            nextgen.append(child2)
        stats.append(average(fitlist))
        print time.strftime("%H:%M:%S"), "Div: ", len(GroupSpecies(gen)), "Ave: ", average(fitlist), "Best: ", fitlist.max()
        #print time.strftime("%H:%M:%S"), "Ave: ", average(fitlist), "Best: ", fitlist.max()
        gen = nextgen # switch local gen variable to point to the next generation

    return gen,stats

# Direct algorithm analysis    
def TestSelection(fitlist):
    sel = []
    sortfit = sort(fitlist) #check if this changes fitlist!
    
    for i in arange(100*Npop):
        sel.append(Selection(sortfit))
    
    hist(sel)

def TestFitness(fitlist):
# answers what percent std(fitlist) is of mean(fitlist). 
# It should be close to 100% for good variability.
    return 100.0*std(fitlist)/mean(fitlist)

# Diversity analysis    
def Diversity(gen):
    div=[]
    for i in range(Npop):
        for j in range(Npop):
            if i != j:
                div.append(CorrDistance(gen[i],gen[j]))
    hist(div,bins=30)
    print "Min: ", min(div), " Max: ", max(div), " Ave: ", mean(div), " Std: ", std(div)

def CorrDiversity(gen):
    div=[]
    for i in range(Npop):
        for j in range(Npop):
            if i != j:
                div.append(Correlation(gen[i],gen[j]))
    hist(div,bins=30)
    print "Min: ", min(div), " Max: ", max(div), " Ave: ", mean(div), " Std: ", std(div)


# Density of states analysis
def DegeneracySpectrum(gen):
    deg=[]
    for i in range(Npop):
        deg.append(LowDOS(gen[i]))
    hist(deg,bins=30)
    print "Min: ", min(deg), " Max: ", max(deg), " Ave: ", mean(deg), " Std: ", std(deg)

def DOS(J):
    tildeJ = FourierTransform(J)
    
    lams = []
    for n1 in range(N1):
        for n2 in range(N2):
            lams.extend(eigvalsh(tildeJ[:,:,n1,n2]))
    
    lams.sort()
    hist(lams)
    
def TestHermetian(h):
    fJ=h.FourierTransform()
    
    total=0.0
    for n1 in range(N1):
        for n2 in range(N2):
            t = abs(fJ[:,:,n1,n2]-fJ[:,:,n1,n2].conj().T).sum()
            total = total + t 
    return total

def TestAllHermetian(gen):
    print "Testing All Hermetian"
    for g in gen:
        print TestHermetian(g)

def Z2Trans(h):
    K=Ham(h)
    
    for a in range(9):
        for b in range(9):
            for m1 in range(N1):
                for m2 in range(N2):
                    K.J[a,b,m1,m2] = h.J[b,a,-m1%N1,-m2%N2]
    return K