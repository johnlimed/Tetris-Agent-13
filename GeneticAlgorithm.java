import java.util.*;
// import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

public class GeneticAlgorithm {
public int populationSize;
public float crossoverRate; 
public static final int NUM_GAMES = 5; // number of games to run to assess fitness of an individual
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

// incomplete
public void train(int generations) {
	assert(generations > 0);
	int bestIndex = 0;
FitnessAssessment bestFitness = assessFitness(population.get(0));

for (int individual = 1; individual < population.size(); individual++) {
	FitnessAssessment fitness = assessFitness(population.get(individual));
	
	if (bestFitness.compareTo(fitness)  > 0) {
		bestFitness = fitness;
		bestIndex = individual;
	}
}


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
		GeneticAlgorithm ga = new GeneticAlgorithm(100);
		
		// Fitness function aka happiness function: aggregate of the 5 heuristics we are using, rank original states by highest happiness
		// Selection: choose 2 parent states at random from a pool of the top 30% fittest in the population to mate
		// Crossover: select a random point to mix and match between 2 parents to get 2 offspring
		// Mutation: select a random point from each off spring to mutate
		// calculate the fitness function of each offspring
		// return the offspring with the highest happiness
		// end of GA

		//***************************************** HELPER FUNCTIONS *************************************************

		// fitness function: to be implmented
		// selection function: to be implemented
		// crossover: to be implemented
		// mutation: to be implemented



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
	}

	
}
