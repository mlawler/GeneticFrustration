module kagome

using Graphs
import Base.*
import Base.+

export symJ,asymJ,nnJ,ddJ,symTerms,asymTerms
#export labeltosite,sitetoindex,indextosite,rtopos
export Ham, GeneticAlgorithm, tograph, Bands
export Initialization,Fitness,Selection,Crossover,Mutation

N1 = 10
N2 = 10
Npop = 100
mutation_rate=0.01 #careful with this: a high mutation rate can saturate the population with mutants!
Ngen=1000
lamfrac=0.08

# Species identification?
corsig(n) = stdev([corr(rand(n),rand(n)) for i in 1:100])
SpeciesCorr=0.1 # A correlation coefficient greater than this is significant
SpeciesPop=20

########## Definition of species ##################
# a species is an object! It has well defined rules for what it can
# do and what its properties are.

# Note: __X is how you declare a private variable. Can we declare private functions too?

# J1-J2-J3-J3p Heisenberg model
symJ= (1.0,0.2,0.01,0.01) # The DNA for a particular J1-J2-J3-J3p model
asymJ=ones(6+6+3+6)
asymJ[7:12]=0.2*ones(6)
asymJ[13:15]=0.01*ones(3)
asymJ[16:end]=0.01*ones(6)
nnJ = [i < 7 ? 1.0 : 0.0 for i in 1:21]
ddJ = ones(6+6+3+6)
ddJ[7:12] = 0.2*ones(6)
ddJ[13:15] = zeros(3)
ddJ[16:end] = 0.2*ones(6)

########## Structures of the Hamiltonian ##############
# These should be moved to a file (eventually)

symTerms=Array{Array{Tuple{Int64,Int64,Int64,Int64},1},1}([])
asymTerms=Array{Array{Tuple{Int64,Int64,Int64,Int64},1},1}([])

# J1 term
push!(symTerms,[(1,2,0,0),(1,3,0,0),(2,1,0,0),(2,3,0,0),(3,1,0,0),(3,2,0,0),
                (2,1,1,0),(3,2,-1,0),(3,2,0,1),(2,3,0,-1),(3,1,1,1),(1,3,-1,-1)])

#=push!(symTerms,[(1,4,0,0),(2,5,0,0),(3,6,0,0),(1,7,0,0),(2,8,0,0),(3,9,0,0),
             (4,1,0,0),(5,2,0,0),(6,3,0,0),(4,7,0,0),(5,8,0,0),(6,9,0,0),
             (7,1,0,0),(8,2,0,0),(9,3,0,0),(7,4,0,0),(8,5,0,0),(9,6,0,0),
             (4,1,1,0),(5,2,1,0),(6,3,1,0),(7,4,-1,0),(2,5,-1,0),(3,6,-1,0),
             (7,4,0,1),(8,5,0,1),(9,6,0,1),(4,7,0,-1),(5,8,0,-1),(6,9,0,-1),
             (7,1,1,1),(8,2,1,1),(9,3,1,1),(1,7,-1,-1),(2,8,-1,-1),(3,9,-1,-1)])=#

# J2 term
push!(symTerms,[(3,1,1,0),(2,3,1,0),(1,3,-1,0),(3,2,-1,0),(1,2,0,1),(3,1,0,1),
                (1,3,0,-1),(2,1,0,-1),(2,1,1,1),(3,2,1,1),(1,2,-1,-1),(2,3,-1,-1)])

#=push!(symTerms,[(7,1,1,0),(8,2,1,0),(9,3,1,0),(4,7,1,0),(5,8,1,0),(6,9,1,0),
             (1,7,-1,0),(2,8,-1,0),(3,9,-1,0),(7,4,-1,0),(8,5,-1,0),(9,6,-1,0),
             (1,4,0,1),(2,5,0,1),(3,6,0,1),(7,1,0,1),(8,2,0,1),(9,3,0,1),
             (1,7,0,-1),(2,8,0,-1),(3,9,0,-1),(4,1,0,-1),(5,2,0,-1),(6,3,0,-1),
             (4,1,1,1),(5,2,1,1),(6,3,1,1),(7,4,1,1),(8,5,1,1),(9,6,1,1),
             (1,4,-1,-1),(2,5,-1,-1),(3,6,-1,-1),(4,7,-1,-1),(5,8,-1,-1),(6,9,-1,-1)])=#
# J3 terms
push!(symTerms,[(3,3,1,0),(3,3,-1,0),(1,1,0,1),(1,1,0,-1),(2,2,1,1),(2,2,-1,-1)])
#=push!(symTerms,[(7,7,1,0),(8,8,1,0),(9,9,1,0),(7,7,-1,0),(8,8,-1,0),(9,9,-1,0),
             (1,1,0,1),(2,2,0,1),(3,3,0,1),(1,1,0,-1),(2,2,0,-1),(3,3,0,-1),
             (4,4,1,1),(5,5,1,1),(6,6,1,1),(4,4,-1,-1),(5,5,-1,-1),(6,6,-1,-1)])=#
# J3p terms
push!(symTerms,[(1,1,1,0),(2,2,1,0),(1,1,-1,0),(2,2,-1,0),(2,2,0,1),(3,3,0,1),
                (2,2,0,-1),(3,3,0,-1),(1,1,1,1),(3,3,1,1),(1,1,-1,-1),(3,3,-1,-1)])

#=push!(symTerms,[(1,1,1,0),(2,2,1,0),(3,3,1,0),(4,4,1,0),(5,5,1,0),(6,6,1,0),
             (1,1,-1,0),(2,2,-1,0),(3,3,-1,0),(4,4,-1,0),(5,5,-1,0),(6,6,-1,0),
             (4,4,0,1),(5,5,0,1),(6,6,0,1),(7,7,0,1),(8,8,0,1),(9,9,0,1),
             (4,4,0,-1),(5,5,0,-1),(6,6,0,-1),(7,7,0,-1),(8,8,0,-1),(9,9,0,-1),
             (1,1,1,1),(2,2,1,1),(3,3,1,1),(7,7,1,1),(8,8,1,1),(9,9,1,1),
             (1,1,-1,-1),(2,2,-1,-1),(3,3,-1,-1),(7,7,-1,-1),(8,8,-1,-1),(9,9,-1,-1)])=#

# J1 terms
push!(asymTerms,[(1,2,0,0),(2,1,0,0)])
push!(asymTerms,[(1,3,0,0),(3,1,0,0)])
push!(asymTerms,[(2,3,0,0),(3,2,0,0)])
push!(asymTerms,[(2,1,1,0),(1,2,-1,0)])
push!(asymTerms,[(3,2,0,1),(2,3,0,-1)])
push!(asymTerms,[(3,1,1,1),(1,3,-1,-1)])

# J2 terms
push!(asymTerms,[(3,1,1,0),(1,3,-1,0)])
push!(asymTerms,[(2,3,1,0),(3,2,-1,0)])
push!(asymTerms,[(1,2,0,1),(2,1,0,-1)])
push!(asymTerms,[(3,1,0,1),(1,3,0,-1)])
push!(asymTerms,[(2,1,1,1),(1,2,-1,-1)])
push!(asymTerms,[(3,2,1,1),(2,3,-1,-1)])

# J3 terms
push!(asymTerms,[(3,3,1,0),(3,3,-1,0)])
push!(asymTerms,[(1,1,0,1),(1,1,0,-1)])
push!(asymTerms,[(2,2,1,1),(2,2,-1,-1)])

# J3p terms
push!(asymTerms,[(1,1,1,0),(1,1,-1,0)])
push!(asymTerms,[(2,2,1,0),(2,2,-1,0)])
push!(asymTerms,[(2,2,0,1),(2,2,0,-1)])
push!(asymTerms,[(3,3,0,1),(3,3,0,-1)])
push!(asymTerms,[(1,1,1,1),(1,1,-1,-1)])
push!(asymTerms,[(3,3,1,1),(3,3,-1,-1)])

function labeltosite(term)
   a,b,m1,m2 = term
   site1=(0, 0, a-1) #(spin component, m1, m2, sublattice index)
   site2=(m1, m2, b-1) #(spin component, m1, m2, sublattice index)
   return site1,site2
end

function sitetoindex(site)
   m1,m2,mu = site
   return mod(m1,N1)+mod(m2,N2)*N1+mu*N1*N2+1
end

function indextosite(i)
   mu = floor(Int,(i-1)/N1/N2)
   m2 = floor(Int,(i-mu*N1*N2-1)/N1)
   m1 = i-mu*N1*N2-m2*N1-1
   return (m1,m2,mu)
end

a1 = (1.0,0.0)
a2 = (-1/2,sqrt(3)/2)
d = [(0.0,0.0),(0.5,0.0),(0.25,sqrt(3)/4)]

*(a::Real, r::Tuple{Float64,Float64}) = (a*r[1],a*r[2])
+(r1::Tuple{Float64,Float64},r2::Tuple{Float64,Float64}) = (r1[1]+r2[1],r1[2]+r2[2])

function rtopos(site)
   m1,m2,mu = site
   return m1*a1 + m2*a2 + d[mu+1]
end

type Ham
  model::Array{Array{Tuple{Int64,Int64,Int64,Int64},1},1}
  genome::Array{Float64,1}

  function Ham(m,g=[])
    if length(g) == 0
      new(m,rand(length(m))-0.5)
    else
      new(m,g)
    end
  end
end

function Base.show(io::IO, h::Ham)
  format(sa) = if sa[1]=='-'
                  return sa[1:min(5,length(sa))]
               else
                  return sa[1:min(4,length(sa))]
               end

  out=format(string(h.genome[1]))
  for a in h.genome[2:end]
     if abs(a) < 0.0001
        out=out*" "*format(string(0.0))
     else
        out=out*" "*format(string(a))
     end
  end

  print(io, out)
end

function FourierTransform(h::Ham,fN1::Int64=N1,fN2::Int64=N2)
# Lesson from profiling this function is that multidimensional arrays have
# slow indexing access. Below, I minimize the number of times I index tildeJ.
  tildeJ = zeros(Complex{Float64},(3,3,fN1,fN2))

  for n1 in 1:fN1
     for n2 in 1:fN2
        temp = zeros(Complex{Float64},(3,3))
        for i in 1:length(h.model)
           for term in h.model[i]
              a,b,m1,m2 = term
              temp[a,b] += h.genome[i]*exp(-1.0im*2.0*pi*(m1*n1/fN1+m2*n2/fN2))
           end
         end
         tildeJ[:,:,n1,n2] = temp
     end
  end
  return tildeJ
end

function LowDOS(h::Ham)
    tildeJ = FourierTransform(h)

    lams = []
    for n1 in 1:N1
        for n2 in 1:N2
            append!(lams,eigvals(tildeJ[:,:,n1,n2]))
        end
    end

    sort!(lams)
    lammin = lams[1]
    lammax = lams[end]
    W = lammax-lammin  # Bandwidth

    dos = 0
    for lam in lams
        if lam-lammin < lamfrac*W #lamfrac is a global variable here
            dos = dos + 1
        end
    end

    return dos #/ (1.0*9*N1*N2)
end

############ Initialization ####################

function Initialization(terms=asymTerms, Np=Npop)
    igen=Array{Ham,1}([])
    for i in 1:Np
        push!(igen,Ham(terms))
    end
    return igen
end

############# Fitness ####################

function Fitness(gen,microHam::Ham)
    return IndividualFitness(gen,microHam)
end

function IndividualFitness(gen,microHam::Ham)
    Np = length(gen)
    fitlist = zeros(Np)
    for i in 1:Np
        dist=Distance(gen[i],microHam)
        lDOS = LowDOS(gen[i])
        fitlist[i] = lDOS/(1.0+dist) # Ensure that dist < 1.0 is close enough for search
    end
    return fitlist
end

function Distance(h1::Ham,h2::Ham)
    AnticorrDist = 2.0 # The distance of a perfectly negative correlation
    corr = cor(h1.genome,h2.genome)
    if corr > 0.0
        return 1/corr-1
    elseif corr < 0.0
        return -AnticorrDist/corr
    else
        return inf
    end
end

function GroupSpecies(gen)
    ungrouped = copy(gen)
    grouped = Array{Array{Ham,1},1}([])

    function FindGroup(ungrouped)
        group=Array{Ham,1}([])
        remaining=Array{Ham,1}([])
        example=pop!(ungrouped)
        push!(group,example)
        while length(ungrouped) > 0
            test=pop!(ungrouped)
            if cor(example.genome,test.genome) > SpeciesCorr
                push!(group,test)
            else
                push!(remaining,test)
            end
        end
        return group, remaining
    end

    while length(ungrouped) > 0
        group,remaining = FindGroup(ungrouped)
        push!(grouped,group)
        ungrouped=remaining
    end

    return grouped
end

################ Selection ##################

function Selection(fitlist)
    return RouletteSelection(fitlist)
end

function RouletteSelection(fitlist)
    flcumsum=cumsum(fitlist)
    p=rand()*flcumsum[end]
    for i in 1:length(fitlist)
        if p<flcumsum[i]
            return i
        end
    end
end

################### Crossover ###################

function Crossover(h1::Ham, h2::Ham)
   i=rand(1:length(h1.genome))
   newgenome = copy(h1.genome)
   newgenome[i:end] = h2.genome[i:end]
   return newgenome
end

################## Mutation ####################

function Mutation(genome)
    return ProportionalFlipMutation(genome)
end

function ProportionalFlipMutation(genome)
    newgenome = copy(genome)
    for i in 1:length(genome)
       if rand() < mutation_rate
          newgenome[i] = genome[i]*4.0*(rand()-0.5)
       end
    end
    return newgenome
end

############## Main genetic alogirthm ##################

function GeneticAlgorithm(gen,microHam,stats=Array{Float64}([]))
    # Recording the fitness of the previous generation seems confusing.
    t=time()
    Np = length(gen)
    fitlist=Fitness(gen,microHam)
    for g in 1:Ngen
        nextgen = Array{Ham,1}([])
        for i in 1:Np
            parent1 = Selection(fitlist)
            parent2 = Selection(fitlist)
            igenome = Crossover(gen[parent1],gen[parent2])
            mgenome = Mutation(igenome)
            push!(nextgen,Ham(gen[parent1].model,mgenome))
        end
        gen = nextgen # switch local gen variable to point to the next generation
        fitlist=Fitness(gen,microHam)

        # Print stats
        push!(stats,mean(fitlist))
        println("$(floor(time()-t)) Div: $(length(GroupSpecies(gen))), Ave: $(mean(fitlist)), Best: $(maximum(fitlist))")
    end
    return gen,stats
end

################## Analyzing results ######################
function tograph(h::Ham)
   # h is the hamiltonian
   # alpha0,beta0 are the selection of alpha, beta in J_{ialpha,jbeta}

   g = simple_graph(0)  # Start with an empty graph

   v = []
   edgewidth = []
   for t in 1:length(h.model)
      for term in h.model[t]
         s1,s2=labeltosite(term)
         (x1,y1,mu1) = s1
         (x2,y2,mu2) = s2
         R =(N1/2,N2/2,0)  # shift origin to middle of system (m1,m2,mu)
         i=sitetoindex((x1+R[1],y1+R[2],mu1+R[3]))
         j=sitetoindex((x2+R[1],y2+R[2],mu2+R[3]))

         clusteri = findfirst(v,i)
         clusterj = findfirst(v,j)

         # I guess the following assumes clusteri != clusterj
         if clusteri == 0
            push!(v,i)
            add_vertex!(g)
            clusteri = findfirst(v,i)
         end

         if clusterj == 0
            push!(v,j)
            add_vertex!(g)
            clusterj = findfirst(v,j)
         end

         add_edge!(g,clusteri,clusterj)
         push!(edgewidth,h.genome[t])
      end
   end

   Nv = length(v)
   x = zeros(Nv)
   y = zeros(Nv)
   for i in 1:Nv
      x[i],y[i]=rtopos(indextosite(v[i]))
   end

   return g,x,y,v,edgewidth
end

function Bands(h::Ham)
   fN1=100
   fN2=100
   lam=zeros((3,fN1,fN2))
   tildeJ = FourierTransform(h,fN1,fN2)
   fshift1(n1) = mod(n1-1+floor(Int,fN1/2),fN1)+1
   fshift2(n2) = mod(n2-1+floor(Int,fN2/2),fN2)+1

   lammax = 0 # needed to compute bandwidth
   for n1 in 1:fN1
      for n2 in 1:fN2
         lams = eigvals(tildeJ[:,:,n1,n2])
         if maximum(lams) > lammax
            lammax=maximum(lams)
         end
         lam[:,fshift1(n1),fshift2(n2)] = sort(lams)
      end
   end
   return (lam-minimum(lam))/(lammax - minimum(lam))
end


end # Module GenFrustration
