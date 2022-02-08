#!/usr/bin/env python3

# 2-D cellular automata simulating tumor growth
# Jessen Havill, Denison University

import turtle
import random
import copy
import numpy
import matplotlib.pyplot as mpl
from heapq import *

COLUMNS = 120
ROWS = 120
SCALE = 10

# Default values

TELOMERE_LENGTH = 50         # tl
EVADE_APOPTOSIS = 10         # e
BASE_MUTATION_RATE = 1000    # m
GENETIC_INSTABILITY = 100    # i
IGNORE_GROWTH_INHIBIT = 10   # g
RANDOM_APOPTOSIS = 1000      # a

# Values to facilitate cancer

# TELOMERE_LENGTH = 35         # tl
# EVADE_APOPTOSIS = 20         # e
# BASE_MUTATION_RATE = 1000    # m
# GENETIC_INSTABILITY = 100    # i
# IGNORE_GROWTH_INHIBIT = 4    # g
# RANDOM_APOPTOSIS = 400       # a

NEIGHBORS = [(-1, -1), (-1, 0), (-1, 1), (1, -1), (1, 0), (1, 1), (0, -1), (0, 1)]

neighbors = { }      # remember previously computed neighborhoods

DRAW = True          # whether to draw CA grid during simulation
TIME = 750           # total simulation time
TRIALS = 1           # number of trials

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
                
def drawSquare(pos, color, tortoise):
    if DRAW:
        (row, column) = pos
        row = ROWS - row - 1
        tortoise.shape(color)
        tortoise.up()
        tortoise.goto(column, row + 1)
        tortoise.stamp()
    
def drawGrid(tortoise):
    if DRAW:
        tortoise.pencolor('gray')
        for row in range(ROWS + 1):
            tortoise.up()
            tortoise.goto(0, row)
            tortoise.down()
            tortoise.goto(COLUMNS, row)
        for column in range(COLUMNS + 1):
            tortoise.up()
            tortoise.goto(column, 0)
            tortoise.down()
            tortoise.goto(column, ROWS)
        
def createSquares(screen, colors):
    square = ((0, 0), (0, SCALE), (SCALE, SCALE), (SCALE, 0))
    for color in colors:
        squareShape = turtle.Shape('compound')
        squareShape.addcomponent(square, color, 'gray')
        screen.register_shape(color, squareShape)
    
def getNeighbors(grid, position):
    global neighbors
    if position in neighbors:
        return neighbors[position]
    neighborly = []
    for n in NEIGHBORS:
        c = position[1] + n[1]
        r = position[0] + n[0]
        if c in range(COLUMNS) and r in range(ROWS):
            neighborly.append((r,c))
    neighbors[position] = neighborly
    return neighborly
    
def emptyNeighbors(grid, position):
    empty = []
    for n in getNeighbors(grid, position):
        if grid[n[0]][n[1]] == None:
            empty.append(n)
    return empty
    
def die(grid, position, tortoise, healthy, mutants):
    cell = grid[position[0]][position[1]]
    if cell.cancerous:
        mutants -= 1
    else:
        healthy -= 1
    del cell
    grid[position[0]][position[1]] = None
    drawSquare(position, 'white', tortoise)
    
    return healthy, mutants
    
def grow(end, tortoise):
    drawGrid(tortoise)
    
    grid = []                    # grid[x][y] = cell at this position or None
    for r in range(ROWS):
        row = [None] * COLUMNS
        grid.append(row)
        
    center = (ROWS // 2, COLUMNS // 2)
    grid[ROWS // 2][COLUMNS // 2] = Cell()     # one healthy cell at center
    drawSquare(center, 'darkgreen', tortoise)
    
    clock = 0               # initialize clock
    pq = [(0, center)]      # insert first mitosis event in priority queue
    
    cancerCounts = numpy.zeros(end)   # cancer cell count at time t
    healthyCounts = numpy.zeros(end)  # healthy cell count at time t
    mutants = 0                       # current cancer cell count
    healthy = 1                       # current healthy cell count
        
    while (len(pq) > 0) and (clock < end):
        (time, position) = heappop(pq)
            
        if time > clock:                             # clock is advancing
            healthyCounts[clock] = healthy           # record cell counts
            cancerCounts[clock] = mutants
        
        clock = time      # update clock with time of new event
        
        cell = grid[position[0]][position[1]]   # get cell at event position
        if cell == None:                        # cell died
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
            healthy, mutants = die(grid, position, tortoise, healthy, mutants)
            continue
                
        # Check for mitosis.
         
        mitosis = True
        empty = emptyNeighbors(grid, position)
        if not cell.selfGrowth():                     # tissue extent check
            if position[0] < 0.05 * ROWS or \
               position[0] > 0.95 * ROWS or \
               position[1] < 0.05 * COLUMNS or \
               position[1] > 0.95 * COLUMNS:
                mitosis = False
        if mitosis:
            if not cell.ignoreAntigrowth():           # ignore growth inhibit check
                if len(empty) == 0:
                    mitosis = False
            elif random.random() > 1 / IGNORE_GROWTH_INHIBIT:
                mitosis = False
                
        # Perform mitosis if warranted.
        
        if mitosis:
        
            # Choose random cell to divide into.
            if len(empty) == 0:       # ignore anti-growth hallmark is on
                daughterPosition = random.choice(getNeighbors(grid, position))
                healthy, mutants = die(grid, daughterPosition, tortoise, healthy, mutants)
            else:
                daughterPosition = random.choice(empty)
                
            cell.telomereLength -= 1           # decrease telomere length
            
            daughter = copy.deepcopy(cell)     # create daughter cell
            
            # Mutate.
            if daughter.geneticInstability():  # increase mutation rate
                daughter.mutationRate = BASE_MUTATION_RATE / GENETIC_INSTABILITY
            if random.random() < 1 / daughter.mutationRate:
                daughter.mutate()
            
            grid[daughterPosition[0]][daughterPosition[1]] = daughter   # insert new cell into grid
            
            if daughter.cancerous:                                      # update cell counts and draw
                drawSquare(daughterPosition, 'yellow', tortoise)
                mutants += 1
            else:
                drawSquare(daughterPosition, 'darkgreen', tortoise)
                healthy += 1
                
            # Schedule mitosis for daughter cell.
            heappush(pq, (clock + random.randrange(5, 11), daughterPosition))
            
        # Schedule mitosis for parent cell.
        heappush(pq, (clock + random.randrange(5, 11), position))
        
    return healthyCounts, cancerCounts
    
def main():
    totalHealthyCounts = numpy.zeros(TIME)       # init total counts for all trials
    totalCancerCounts = numpy.zeros(TIME)
    
    for t in range(TRIALS):                      # run TRIALS trials
        if DRAW:
            tortoise = turtle.Turtle()
            screen = tortoise.getscreen()
            screen.setup(COLUMNS * SCALE + 20, ROWS * SCALE + 20)
            screen.setworldcoordinates(0, 0, COLUMNS, ROWS)
            screen.tracer(500)
            tortoise.hideturtle()
            createSquares(screen, ['darkgreen', 'yellow', 'white'])
        else:
            tortoise = None
    
        healthyCounts, cancerCounts = grow(TIME, tortoise)
        totalHealthyCounts += healthyCounts
        totalCancerCounts += cancerCounts
    
        if DRAW:
            screen.update()
            screen.exitonclick()
    
    # Display plot of cell counts over time.
    
    totalHealthyCounts /= TRIALS
    totalCancerCounts /= TRIALS
    mpl.plot(range(TIME), totalHealthyCounts, label = 'Healthy')
    mpl.plot(range(TIME), totalCancerCounts, label = 'Cancerous')
    mpl.legend(loc = 'upper left')
    mpl.xlabel('Time')
    mpl.ylabel('Number of cells')
    mpl.show()
    
main()
