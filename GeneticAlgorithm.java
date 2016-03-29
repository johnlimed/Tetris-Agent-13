import java.util.*;
// import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

public class GeneticAlgorithm {
	private static final boolean isLogging = true;
	private static final Logger logger = Logger.getLogger("GeneticAlgorithm");
	public int populationSize;
	public static float CROSSOVER_RATE = 0.5f;
	public static final int NUM_GAMES = 5; // number of games to run to assess fitness of an individual
	public static final int TOURNAMENT_SIZE = 2; // 2's the most common setting. 1 is random selection, higher values causes higher selection pressure
	ArrayList<ArrayList<FeatureWeightPair>> population;
	PlayerSkeleton player;
	private Random rng;

	public GeneticAlgorithm(int populationSize) {
		rng = new Random();
		player = new PlayerSkeleton();
		population = new ArrayList<ArrayList<FeatureWeightPair>>(populationSize);

		for (int i = 0; i < populationSize; i++) {
			population.add(generateRandomIndividual());
		}
		log("Random individuals generated."); 
		logPopulation();
	}

	private void logPopulation() {
		if (!isLogging)
			return;
		
String str = "population of " + population.size() + " individuals are\r\n";

for (int i=0; i<population.size(); i++) {
		str += i + ": " + getIndividualAsStr(i);
		str += "\r\n";
}
		
log(str);
	}
	
	private String getIndividualAsStr(ArrayList<FeatureWeightPair> individual) {
		String str = "";
		
		for (FeatureWeightPair f : individual) {
			str += f.weight + ", ";
		}
		
		return str;
	}
	
	private String getIndividualAsStr(int i) {
		return getIndividualAsStr(population.get(i));
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
		log("training for " + generations + " generations:");
		
		for (int generation = 0; generation < generations; generation++) {
			System.out.println("currently on generation " + generation);
			log("currently on generation " + generation);
			population = reproduce();
			log("population after reproduction:");
			logPopulation();
		}

		return findBestIndividual();
	}

	// reproduces children
	private ArrayList<ArrayList<FeatureWeightPair>> reproduce() {
		int iterations = population.size()/2;
		ArrayList<ArrayList<FeatureWeightPair>> children = new ArrayList<ArrayList<FeatureWeightPair>>(population.size()); // the next generation

		for (int i = 0; i<iterations; i++) {
			// find 2 parents to mate
			ArrayList<FeatureWeightPair> child1 = deepCopyIndividual(tournamentSelection(TOURNAMENT_SIZE));
			ArrayList<FeatureWeightPair> child2 = deepCopyIndividual(tournamentSelection(TOURNAMENT_SIZE));
			
			if (isLogging)
			log("mating " + getIndividualAsStr(child1) + " with " + getIndividualAsStr(child2));
			
			uniformCrossover(child1, child2);
			
			if (isLogging)
			log("After crossover: child1 = " + getIndividualAsStr(child1) + ", child2 = " + getIndividualAsStr(child2));

			mutate(child1);
			mutate(child2);
			
			if (isLogging)
				log("After mutation: child1 = " + getIndividualAsStr(child1) + ", child2 = " + getIndividualAsStr(child2));
		
			children.add(child1);
			children.add(child2);
		}

		return children;
	}

	// does a deep copy of an individual's weights
	ArrayList<FeatureWeightPair> deepCopyIndividual(ArrayList<FeatureWeightPair> individual) {
		ArrayList<FeatureWeightPair> copy = new ArrayList<FeatureWeightPair>(individual.size());

		for (int i = 0; i < individual.size(); i++)
			copy.add(new FeatureWeightPair(individual.get(i).feature, individual.get(i).weight));

		return copy;
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
		int lowest = player.playGame(false);
		int highest = lowest, total = lowest;

		for (int game = 1; game < NUM_GAMES; game++) {
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

	// crosses over 2 individuals using uniform crossover
	private void uniformCrossover(ArrayList<FeatureWeightPair> x, ArrayList<FeatureWeightPair> y) {
		for (int i = 0; i < x.size(); i++) {
			if (CROSSOVER_RATE >= randomFloat(0.0f, 1.0f)) {
				// swap the genes on these 2 vectors
				float temp = x.get(i).weight;
				x.get(i).weight = y.get(i).weight;
				y.get(i).weight = temp;
			}
		}
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

	public static void loggerInit() {
		try {
			FileHandler handler = new FileHandler("GeneticAlgorithm.txt");
			handler.setFormatter(new SimpleFormatter());
			LogManager.getLogManager().reset();
			logger.addHandler(handler);
			logger.log(Level.INFO, "logger initialized\n");
		} catch (Exception e) {
			/* error opening the log file - just get rid of logging so it won't 
			 * print to the console while the user is running the program */
			LogManager.getLogManager().reset();
		}
	}

	private static void log(String msg) {
		if (isLogging)
		logger.log(Level.INFO, msg);
	}
	
	public static void main(String[] args) {
		loggerInit();
	log("test");  	
		GeneticAlgorithm ga = new GeneticAlgorithm(6); // population size
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
