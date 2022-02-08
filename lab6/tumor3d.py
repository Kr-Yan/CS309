#!/Users/havill/opt/anaconda3/bin/python3
####!/usr/bin/env python3

# 3-D cellular automata simulating tumor growth
# Jessen Havill, Denison University

import random
import copy
import numpy
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import Axes3D
from heapq import *

COLUMNS = 50
ROWS = 50
DEPTH = 50

# Default values

# TELOMERE_LENGTH = 100        # tl
# EVADE_APOPTOSIS = 10         # e
# BASE_MUTATION_RATE = 100000  # m
# GENETIC_INSTABILITY = 100    # i
# IGNORE_GROWTH_INHIBIT = 10   # g
# RANDOM_APOPTOSIS = 1000      # a

# Values to facilitate cancer

TELOMERE_LENGTH = 35         # tl
EVADE_APOPTOSIS = 20         # e
BASE_MUTATION_RATE = 100000  # m
GENETIC_INSTABILITY = 100    # i
IGNORE_GROWTH_INHIBIT = 4    # g
RANDOM_APOPTOSIS = 400       # a

NEIGHBORS = [(-1, -1, -1), (-1, -1, 0), (-1, -1, 1), (-1, 0, -1), (-1, 0, 0),
             (-1, 0, 1), (-1, 1, -1), (-1, 1, 0), (-1, 1, 1), (0, -1, -1), (0, -1, 0),
             (0, -1, 1), (0, 0, -1), (0, 0, 1), (0, 1, -1), (0, 1, 0), (0, 1, 1),
             (1, -1, -1), (1, -1, 0), (1, -1, 1), (1, 0, -1), (1, 0, 0), (1, 0, 1),
             (1, 1, -1), (1, 1, 0), (1, 1, 1)]

neighbors = { }      # remember previously computed neighborhoods

TIME = 520                                # total simulation time
DISPLAY_TIMES = (100, 200, 300, 400, 420, 440, 460, 480, 500)  # times to display growth
DISPLAY_FRACTION = 0.1                    # percent quality of growth displays
TRIALS = 1                                # number of trials

NUM_HALLMARKS = 5
SG = 0            # self-growth
IGI = 1           # ignore antigrowth signals
EA = 2            # evade apoptosis
EI = 3            # effective immortality
GI = 4            # genetic instability

# AG -- ability for angiogenesis -- NOT IMPLEMENTED
# MT -- metastasis -- NOT IMPLEMENTED

class Cell:
    def __init__(self):        
        self.hallmarks = [False] * NUM_HALLMARKS
        self.mutations = 0
        self.cancerous = False                    # True if any hallmark is present
        self.telomereLength = TELOMERE_LENGTH
        self.mutationRate = BASE_MUTATION_RATE
        
    def selfGrowth(self):
        return self.hallmarks[SG]
        
    def ignoreAntigrowth(self):
        return self.hallmarks[IGI]
        
    def evadeApoptosis(self):
        return self.hallmarks[EA]
        
    def effectiveImmortality(self):
        return self.hallmarks[EI]
        
    def geneticInstability(self):
        return self.hallmarks[GI]
        
    def mutate(self):
        if self.mutations < NUM_HALLMARKS:
            j = random.randrange(NUM_HALLMARKS)
            while self.hallmarks[j]:
                j = random.randrange(NUM_HALLMARKS)
            self.hallmarks[j] = True
            self.mutations += 1
            self.cancerous = True
            
def getNeighbors(grid, position):
    global neighbors
    if position in neighbors:
        return neighbors[position]
    neighborly = []
    for n in NEIGHBORS:
        nPosition = (position[0] + n[0], position[1] + n[1], position[2] + n[2])
        if nPosition[0] in range(-ROWS//2, ROWS//2) and nPosition[1] in range(-COLUMNS//2, COLUMNS//2) and nPosition[2] in range(-DEPTH//2, DEPTH//2):
            neighborly.append(nPosition)
    neighbors[position] = neighborly
    return neighborly
    
def emptyNeighbors(grid, position):
    empty = []
    for n in getNeighbors(grid, position):
        if n not in grid:
            empty.append(n)
    return empty
        
def die(grid, position, healthy, mutants):
    cell = grid[position]
    if cell.cancerous:
        mutants -= 1
    else:
        healthy -= 1
    del cell
    del grid[position]
    
    return healthy, mutants
    
def displayGrowth(grid, time):
    fig = mpl.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.set_xlim(-COLUMNS/2, COLUMNS/2)
    ax.set_ylim(-ROWS/2, ROWS/2)
    ax.set_zlim(-DEPTH/2, DEPTH/2)
    ax.set_title('t = ' + str(time))
    mpl.axis('off')
    
    x = []
    y = []
    z = []
    colors = []
    for pos in grid:
        if random.random() < DISPLAY_FRACTION:
            x.append(pos[0])
            y.append(pos[1])
            z.append(pos[2])
            if grid[pos].cancerous:
                colors.append('yellow')
            else:
                colors.append('red')
    
    ax.scatter(x, y, z, c = colors, s = 15)
    mpl.show()
    
def grow(end):
    grid = { }              # grid[(x,y,z)] = cell at this position
        
    center = (0, 0, 0)
    grid[center] = Cell()   # one healthy cell at center
    
    clock = 0               # initialize clock
    pq = [(0, center)]      # insert first mitosis event in priority queue
    
    cancerCounts = numpy.zeros(end)   # cancer cell count at time t
    healthyCounts = numpy.zeros(end)  # healthy cell count at time t
    mutants = 0                       # current cancer cell count
    healthy = 1                       # current healthy cell count
        
    while (len(pq) > 0) and (clock < end):
        (time, position) = heappop(pq)
        
        if time > clock:                          # clock is advancing
            healthyCounts[clock] = healthy        # record cell counts
            cancerCounts[clock] = mutants
            if clock in DISPLAY_TIMES:            # display growth 
                displayGrowth(grid, clock)
            
        clock = time      # update clock with time of new event
        
        cell = grid.get(position, None)     # get cell at event position
        if cell == None:                    # cell died
            continue                        
        
        # Check for cell death.
        
        death = random.random() < 1 / RANDOM_APOPTOSIS              # random apoptosis
        if not cell.evadeApoptosis():
            if random.random() < cell.mutations / EVADE_APOPTOSIS:  # genetic damage
                death = True
        if not cell.effectiveImmortality():
            if cell.telomereLength <= 0:                            # telomere length 0
                death = True
        if death:
            healthy, mutants = die(grid, position, healthy, mutants)
            continue

        # Check for mitosis.

        mitosis = True
        empty = emptyNeighbors(grid, position)
        if not cell.selfGrowth():                           # tissue extent check
            if position[0] < 0.95 * -ROWS/2 or \
               position[0] > 0.95 * ROWS/2 or \
               position[1] < 0.95 * -COLUMNS/2 or \
               position[1] > 0.95 * COLUMNS/2 or \
               position[2] < 0.95 * -DEPTH/2 or \
               position[2] > 0.95 * DEPTH/2:
                mitosis = False
        if mitosis:
            if not cell.ignoreAntigrowth():                 # ignore growth inhibit check
                if len(empty) == 0:
                    mitosis = False
            elif random.random() > 1 / IGNORE_GROWTH_INHIBIT:
                mitosis = False

        # Perform mitosis if warranted.

        if mitosis:

            # Choose random cell to divide into.
            if len(empty) == 0:       # ignore anti-growth hallmark is on
                daughterPosition = random.choice(getNeighbors(grid, position))
                healthy, mutants = die(grid, daughterPosition, healthy, mutants)
            else:
                daughterPosition = random.choice(empty)
                
            cell.telomereLength -= 1           # decrease telomere length
            
            daughter = copy.deepcopy(cell)     # create daughter cell
            
            # Mutate.
            if daughter.geneticInstability():  # increase mutation rate
                daughter.mutationRate = BASE_MUTATION_RATE / GENETIC_INSTABILITY
            if random.random() < 1 / daughter.mutationRate:
                daughter.mutate()
            
            grid[daughterPosition] = daughter  # insert new cell into grid
            
            if daughter.cancerous:             # update cell counts
                mutants += 1
            else:
                healthy += 1
                
            # Schedule mitosis for daughter cell.
            heappush(pq, (clock + random.randrange(5, 11), daughterPosition))
            
        # Schedule mitosis for parent cell.
        heappush(pq, (clock + random.randrange(5, 11), position))
    
    displayGrowth(grid, clock)
    
    return healthyCounts, cancerCounts, grid
    
def main():    
    totalHealthyCounts = numpy.zeros(TIME)        # init total counts for all trials
    totalCancerCounts = numpy.zeros(TIME)
    
    for t in range(TRIALS):                       # run TRIALS trials
        healthyCounts, cancerCounts, grid = grow(TIME)
        totalHealthyCounts += healthyCounts
        totalCancerCounts += cancerCounts
    
    displayGrowth(grid, TIME)                    # display final growth
        
    # Display plot of cell counts over time.
    
    totalHealthyCounts /= TRIALS
    totalCancerCounts /= TRIALS
    fig = mpl.figure()
    mpl.plot(range(TIME), totalHealthyCounts, label = 'Healthy')
    mpl.plot(range(TIME), totalCancerCounts, label = 'Cancerous')
    mpl.legend(loc = 'upper left')
    mpl.xlabel('Time')
    mpl.ylabel('Cell counts')
    mpl.show()
        
main()
