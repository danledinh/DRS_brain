
import numpy as np
import pandas as pd
from sklearn.svm import LinearSVR
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
import random

# --- User configuration ---
INPUT_FEATURES = 'input_features.tsv'  # Tab-separated file, rows: samples, columns: features
INPUT_TARGETS = 'input_targets.tsv'    # Tab-separated file, rows: samples, columns: target variables
TARGET_COLUMN = 'target'               # Name of the target column in INPUT_TARGETS
N_GENERATIONS = 20
POP_SIZE = 30
N_FEATURES = 10  # Number of features to select in each individual
RANDOM_SEED = 42

random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# --- Data loading ---
X = pd.read_csv(INPUT_FEATURES, sep='\t', index_col=0)
y = pd.read_csv(INPUT_TARGETS, sep='\t', index_col=0)[TARGET_COLUMN]
X, y = X.align(y, join='inner', axis=0)
feature_names = list(X.columns)

# --- Fitness function ---
def evaluate_fitness(feature_subset):
    if len(feature_subset) == 0:
        return float('inf')
    X_sub = X[list(feature_subset)]
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_sub)
    svr = LinearSVR(random_state=RANDOM_SEED, max_iter=10000)
    svr.fit(X_scaled, y)
    y_pred = svr.predict(X_scaled)
    mse = mean_squared_error(y, y_pred)
    return mse

# --- Genetic algorithm core ---
def random_individual():
    return tuple(sorted(random.sample(feature_names, N_FEATURES)))

def crossover(parent1, parent2):
    cut = random.randint(1, N_FEATURES-1)
    child = tuple(sorted(set(parent1[:cut] + parent2[cut:])))
    # If child is too small, fill randomly
    while len(child) < N_FEATURES:
        f = random.choice(feature_names)
        if f not in child:
            child = tuple(sorted(child + (f,)))
    # If too large, trim
    child = tuple(sorted(child[:N_FEATURES]))
    return child

def mutate(individual, mutation_rate=0.2):
    ind = list(individual)
    for i in range(N_FEATURES):
        if random.random() < mutation_rate:
            new_f = random.choice([f for f in feature_names if f not in ind])
            ind[i] = new_f
    return tuple(sorted(set(ind)))

# --- Main GA loop ---
population = [random_individual() for _ in range(POP_SIZE)]
fitness_scores = [evaluate_fitness(ind) for ind in population]

for gen in range(N_GENERATIONS):
    # Selection: tournament
    selected = []
    for _ in range(POP_SIZE):
        i, j = random.sample(range(POP_SIZE), 2)
        winner = population[i] if fitness_scores[i] < fitness_scores[j] else population[j]
        selected.append(winner)
    # Crossover and mutation
    next_population = []
    for i in range(0, POP_SIZE, 2):
        parent1, parent2 = selected[i], selected[(i+1)%POP_SIZE]
        child1 = mutate(crossover(parent1, parent2))
        child2 = mutate(crossover(parent2, parent1))
        next_population.extend([child1, child2])
    population = next_population[:POP_SIZE]
    fitness_scores = [evaluate_fitness(ind) for ind in population]
    best_idx = np.argmin(fitness_scores)
    print(f"Generation {gen+1}: Best MSE = {fitness_scores[best_idx]:.4f}")

# --- Output best solution ---
best_idx = np.argmin(fitness_scores)
best_features = population[best_idx]
print("Best feature subset:", best_features)
print("Best MSE:", fitness_scores[best_idx])
