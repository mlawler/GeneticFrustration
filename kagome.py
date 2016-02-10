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

class Ham:
    def __init__(self, oldHam=None, microHam=None, J1=0.0, J2=0.0, J3=0.0, J3p=0.0):
        if oldHam:
            self.J=oldHam.J.copy()
        elif microHam:
            self.J = self.__NonZeroRandomGenerate(microHam.J)
        elif J1 or J2 or J3 or J3p:
            self.J = self.__Microscopic(J1,J2,J3,J3p)
        else:
            self.J = self.__RandomGenerate()
        
    def __RandomGenerate(self):
        J = rand(9,9,N1,N2)-0.5*ones((9,9,N1,N2))
        symJ = zeros(J.shape)
    
        # symmetrize to ensure a Hermetian matrix
        for a in range(9):
            for b in range(9):
                for n1 in range(N1):
                    for n2 in range(N2):
                        symJ[a,b,n1,n2]= (J[a,b,n1,n2]+J[b,a,-n1%N1,-n2%N2])/2.0
        
        return symJ

    def __NonZeroRandomGenerate(self,microJ):
        K=microJ.copy()
        symK = zeros(K.shape)
    
        for a in range(9):
            for b in range(9):
                for n1 in range(N1):
                    for n2 in range(N2):
                        if microJ[a,b,n1,n2] != 0.0:
                            K[a,b,n1,n2]=rand() - 0.5

        # symmetrize to ensure a Hermetian matrix
        for a in range(9):
            for b in range(9):
                for n1 in range(N1):
                    for n2 in range(N2):
                        symK[a,b,n1,n2]= (K[a,b,n1,n2]+K[b,a,-n1%N1,-n2%N2])/2.0

        return symK         

    def __Microscopic(self,J1,J2,J3,J3p): 
        # microscipic model with parameters J1, J2, J3 (cross hex) and J3'
        microJ = zeros((9,9,N1,N2))

        def enlarge(Jmat):
            lJmat = zeros((9,9))
            for i in range(3):
                for j in range(3):
                    lJmat[3*i:3*i+3,3*j:3*j+3] = Jmat[i,j]*eye(3)
            return lJmat
    
        Jmat00 = array([[0,J1,J1],[J1,0,J1],[J1,J1,0]])
        Jmat10 = array([[J3p,0,0],[J1,J3p,J2],[J2,0,J3]])
        Jmatm10 = transpose(Jmat10)
        Jmat01 = array([[J3,J2,0],[0,J3p,0],[J2,J1,J3p]])
        Jmat0m1 = transpose(Jmat01)
        Jmat11 = array([[J3p,0,0],[J2,J3,0],[J1,J2,J3p]])
        Jmatm1m1 = transpose(Jmat11)
    
        microJ[:,:,0,0] = enlarge(Jmat00)
        microJ[:,:,1,0] = enlarge(Jmat10)
        microJ[:,:,-1,0] = enlarge(Jmatm10)
        microJ[:,:,0,1] = enlarge(Jmat01)
        microJ[:,:,0,-1] = enlarge(Jmat0m1)
        microJ[:,:,1,1] = enlarge(Jmat11)
        microJ[:,:,-1,-1] = enlarge(Jmatm1m1)
        
        return microJ

    def FourierTransform(self):
        tildeJ = zeros(self.J.shape,dtype=complex)
    
        for a in range(9):
            for b in range(9):
                tildeJ[a,b] = fft2(self.J[a,b])
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