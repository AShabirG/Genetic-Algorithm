import numpy as np
import matplotlib.pyplot as plt
import random
import bisect
random.seed(4)
n = 10  # Initial number of alleles in gene


# Split genotype to genes
def geno_to_gene(a):
    half = len(a) // 2
    return a[:half], a[half:]


# Function to check if optimal fitness has been achieved in the population
def fitness_check(demes, n):
    for i in range(len(demes)):
        out = optimal_genotype_check(demes[i], n)
        if out:
            return True
    return False


def optimal_genotype_check(deme, n):
    for row in deme:
        geneL, geneR = geno_to_gene(row)
        i = np.count_nonzero(geneL)
        j = np.count_nonzero(geneR)
        if i == n and j == n:
            return True
    return False


# Finds maximum fitness within population and returns index as well as a cumulative sum of fitness in the deme
def loc_of_max_fitness(deme):
    max = 0
    fitnesses = []
    for row in deme:
        geneL, geneR = geno_to_gene(row)
        i = np.count_nonzero(geneL)
        j = np.count_nonzero(geneR)
        fitness = R[i][j] * (2 ** i + 2 ** j)
        fitnesses.append(fitness)
        if fitness > max:
            max = fitness
            max_row = row
    sum_of_fitness = sum(fitnesses)
    fitnesses = np.cumsum(np.true_divide(fitnesses, sum_of_fitness / 100))
    return max_row, fitnesses


# Mutate at a rate of 1 mutation per genotype
def mutation(genotype, n):
    for i in range(len(genotype)):
        if random.randint(1, 2 * n) == 6:  # 1/2n chance of mutation
            if genotype[i] == 1:
                genotype[i] = 0
            else:
                genotype[i] = 1
    return genotype


# One-point crossover implementation
def crossover(parent_1, parent_2):
    cross = random.randint(1, parent_1.shape[0] - 1)
    out = np.concatenate((parent_1[:cross], parent_2[cross:]))
    return out


# Generate noise in the landscape
def generate_R():
    global R, i, j
    R = np.zeros((n + 1, n + 1))
    for i in range(len(R)):
        for j in range(len(R[i])):
            number = 1.5 - random.uniform(0.5, 1)
            R[i][j] = number


# Store n values with the corresponding number of iterations required to solve
n_value = []
iterations_to_solve = []

# For all n values less than 90 in increments of 10
while n < 90:
    optimal = False
    generate_R()  # Generate noise in landscape
    demes = np.random.randint(2, size=(20, 20, 2 * n))  # Initialise 20 demes each containing 20 genotypes of length 2n
    iterations = 0

    while not optimal:
        # Generational reproduction so empty new generation created
        new_gen = np.empty_like(demes)

        # For every deme find the index of maximum fitness to copy over and additionally the cumulative fitness
        fitnesses = []
        for i in range(len(demes)):
            max_row, fit_list = loc_of_max_fitness(demes[i])
            fitnesses.append(fit_list)
            new_gen[i][0] = max_row

        j = 0
        # choose which deme to migrate genotypes between
        chosen_deme = list(range(0, 20))
        random.shuffle(chosen_deme)
        while j < 20:
            # Choose which demes exchange genotypes
            deme_a = random.choice(chosen_deme)
            chosen_deme.remove(deme_a)
            deme_b = random.choice(chosen_deme)
            chosen_deme.remove(deme_b)

            # Randomly choose a value within the cumulative fitness to determine genotype to be migrated
            fitness_a = fitnesses[deme_a][-1] * np.random.rand()
            fitness_b = fitnesses[deme_b][-1] * np.random.rand()

            # Find where the random fitness value falls into the cumulative range
            row_a = bisect.bisect_left(fitnesses[deme_a], fitness_a)
            row_b = bisect.bisect_left(fitnesses[deme_b], fitness_b)

            # Migrate the genotype and mutate it
            new_gen[deme_a][1] = mutation(np.copy(demes[deme_b][row_a]), n)
            new_gen[deme_b][1] = mutation(np.copy(demes[deme_a][row_b]), n)

            j += 2

        for l in range(len(demes)):  # number of demes
            for i in range(2, 20):  # number of blank rows
                # Fitness proportionate selection of remainder of genotypes
                k = fitnesses[l][-1] * np.random.rand()
                m = fitnesses[l][-1] * np.random.rand()
                parent_1 = bisect.bisect_left(fitnesses[l], k)
                parent_2 = bisect.bisect_left(fitnesses[l], m)
                new_gen[l][i] = mutation(crossover(np.copy(demes[l][parent_1]), np.copy(demes[l][parent_2])), n)

        # Check if maximum fitness is in the population
        optimal = fitness_check(new_gen, n)

        # Break condition if taking too many iterations and reattempt the same n value
        if iterations >= 2000:
            print(f"Solution not found in 2000 iterations for n = {n}. Restarting.")
            n -= 10
            break

        iterations += 1
        demes = new_gen

    # if it successfully solved
    if iterations != 2000:
        n_value.append(n)
        iterations_to_solve.append(iterations)
        print(n)
        print(iterations)
    n += 10


print("n=", n_value)
print("iteration=", iterations_to_solve)
plt.plot(n_value, iterations_to_solve)
plt.xlabel("n")
plt.ylabel("iterations")
plt.title("Number of iterations to get optimal without crossover")
plt.show()
