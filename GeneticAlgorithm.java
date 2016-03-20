import java.util.*;
// import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

public class GeneticAlgorithm {
public int populationSize;
public float crossoverRate; 
public static final int NUM_GAMES = 5; // number of games to run to assess fitness of an individual
public static final int TOURNAMENT_SIZE = 2; // 2's the most common setting. 1 is random selection, higher values causes higher selection pressure
ArrayList<ArrayList<FeatureWeightPair>> population;
PlayerSkeleton player;
private Random rng;

public GeneticAlgorithm(int populationSize) {
	rng = new Random();
	player = new PlayerSkeleton();
	population = new ArrayList<ArrayList<FeatureWeightPair>>(populationSize);
	crossoverRate = 1.0f / populationSize;
	
	for (int i = 0; i < populationSize; i++) {
		population.add(generateRandomIndividual());
		/* 
		System.out.print("individual " + i + ":");
		for (FeatureWeightPair f : population.get(i)) {
			System.out.print(f.weight + ", ");
		}
		System.out.print("\n");
		*/
	}
}

private ArrayList<FeatureWeightPair> generateRandomIndividual() {
	ArrayList<FeatureWeightPair> individual = new ArrayList<FeatureWeightPair>();
	// all the feature functions we're using so far contribute negatively to happiness and so should be minimized, 
// 	hence their weights should be negative
	// for example, the presence of holes should decrease happiness
	// not sure if -1.0 is a good lower bound for the initial population
	 	
	individual.add(new FeatureWeightPair(new PlayerSkeleton.AggHeight(), randomFloat(-1.0f, 0.0f)));
	individual.add(new FeatureWeightPair(new PlayerSkeleton.Bumpiness(), randomFloat(-1.0f, 0.0f)));
	individual.add(new FeatureWeightPair(new PlayerSkeleton.MaxHeight(), randomFloat(-1.0f, 0.0f)));
	individual.add(new FeatureWeightPair(new PlayerSkeleton.NumHoles(), randomFloat(-1.0f, 0.0f)));

	return individual;
}

// returns a new float in the range [min, max)
private float randomFloat(float minInclusive, float maxExclusive) {
 return (float) ThreadLocalRandom.current().nextDouble(minInclusive, maxExclusive);
}

//returns a new int in the range [min, max)
// to have both endpoints included, pass a value of max+1 to this function
private int randomInt(int minInclusive, int maxExclusive) {
return  ThreadLocalRandom.current().nextInt(minInclusive, maxExclusive);
}

// run the genetic algorithm for a specified number of generations
// returns fitness information for the best individual in the last generation
public FitnessAssessment trainFor(int generations) {
	assert(generations > 0);
	for (int generation = 0; generation < generations; generation++) {
population = reproduce();
	}
	
	return findBestIndividual();
}

	// reproduces children
	private ArrayList<ArrayList<FeatureWeightPair>> reproduce() {
		int iterations = population.size()/2;
		ArrayList<ArrayList<FeatureWeightPair>> children = new ArrayList<ArrayList<FeatureWeightPair>>(population.size()); // the next generation

		for (int i = 0; i<iterations; i++) {
			ArrayList<FeatureWeightPair> parent1 = tournamentSelection(TOURNAMENT_SIZE);
			ArrayList<FeatureWeightPair> parent2 = tournamentSelection(TOURNAMENT_SIZE);
			// crossover (mate) these parents
			// mutate children
			// add the 2 mutated children to nextGeneration
		}

		return children;
	}
	
	// returns fitness information on the best individual
private FitnessAssessment findBestIndividual() {
	int bestIndex = 0;
FitnessAssessment bestFitness = assessFitness(population.get(0));

for (int individual = 1; individual < population.size(); individual++) {
	FitnessAssessment fitness = assessFitness(population.get(individual));
	
	if (bestFitness.compareTo(fitness)  > 0) {
		bestFitness = fitness;
		bestIndex = individual;
	}
}

	return bestFitness;
}

// returns information about the lowest, average and highest score on an individual after playing NUM_GAMES games, each game with random piece sequences
private FitnessAssessment assessFitness(ArrayList<FeatureWeightPair> individual) {
	player.setFeatureWeightPairs(individual);
	int lowest = 0, highest = 0, total = 0;
	
	for (int game = 0; game < NUM_GAMES; game++) {
		int score = player.playGame(false);
		total += score;
		
		lowest = Math.min(lowest, score);
		highest = Math.max(highest, score);
	}
	
	float average = total * 1.0f / NUM_GAMES;
	return new FitnessAssessment(individual, lowest, average, highest);
}

private ArrayList<FeatureWeightPair> tournamentSelection(int tournamentSize) {
	int best = randomInt(0, population.size());
	FitnessAssessment bestFitness = assessFitness(population.get(best));
	
	for (int i = 2; i <= tournamentSize; i++) {
		int next = randomInt(0, population.size());
		
		// ensure that the next individual selected is different from the existing best
		while (best == next)
			next = randomInt(0, population.size());
		
		FitnessAssessment fitness = assessFitness(population.get(next));
		
		if (bestFitness.compareTo(fitness)  > 0) {
			bestFitness = fitness;
			best = next;
		}
	}
	
	return population.get(best);
}

// crosses over 2 individuals, not implemented 
private void crossover(ArrayList<FeatureWeightPair> x, ArrayList<FeatureWeightPair> y) {

}

// mutates an individual using the Gaussian Convolution algorithm
private void mutate(ArrayList<FeatureWeightPair> individual) {
for (int i = 0; i < individual.size(); i++) {
	float n = 0.0f;
	
	do {
		n  = (float) rng.nextGaussian(); // this uses standard dev of 1.0, not sure if this needs to be changed later
	} while (individual.get(i).weight + n >= 0.0f);
		
		individual.get(i).weight += n;
}

	}

	public static void main(String[] args) {
		GeneticAlgorithm ga = new GeneticAlgorithm(100); // population size
FitnessAssessment result =ga.trainFor(2); // number of generations to train for 		
System.out.println("Training complete. The best individual is ");
System.out.println(result);
				
	}

	// stores information about the fitness of an individual
	public static class FitnessAssessment implements Comparable<FitnessAssessment> {
		public ArrayList<FeatureWeightPair> individual;
		public int lowest, highest; // scores
		public float average; // average score
		
		public FitnessAssessment(ArrayList<FeatureWeightPair> individual, int l, float a, int h) {
			this.individual = individual;
			lowest = l;
			average = a;
			highest = h;
		}
		
		public int compareTo(FitnessAssessment other) {
		if (other.average < average)
			return -1;
		
		if (other.average > average)
			return 1;
		
		// tiebreak using lowest
		if (other.lowest < lowest)
			return -1;
		
		if (other.lowest > lowest)
			return 1;
		
		// then highest
		if (other.highest < highest)
			return -1;
		
		if (other.highest > highest)
			return 1;
		
		return 0;
		}
		
		@Override
		public String toString() {
			String str = "lowest = " + lowest + " average = " + average + " highest = " + highest + "\n";
			str += "weights: ";
			
			for (FeatureWeightPair f : individual)
				str += f.weight + ", ";
			
			return str;
		}
	}

	
}
