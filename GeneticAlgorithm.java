import java.lang.reflect.Array;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
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
	public int populationSize = 0;
	public float crossoverRate = 0.0f;
	public int numGames = 0; // number of games to run to assess fitness of an individual
	public float mutationSigma = 0.0f; // standard deviation for mutation
	public int tournamentSize = 2; // 2's the most common setting. 1 is random selection, higher values causes higher selection pressure
	public int numElites= 4; // number of numElites to keep
	ArrayList<ArrayList<FeatureWeightPair>> population;
	ArrayList<FitnessAssessment> fitnessResults;
	PlayerSkeleton player;
	private Random rng;
	private ExecutorService service = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

	public GeneticAlgorithm(float crossover, int elites, int games, float mutationSigma, int populationSize, int tournamentSize) {
		this.crossoverRate = crossover;
		this.numElites = elites;
		this.numGames = games;
		this.mutationSigma = mutationSigma;
		this.populationSize = populationSize;
		this.tournamentSize = tournamentSize;

		assert(numElites >= 0 && numElites <= populationSize);
		assert ((populationSize - numElites) % 2 == 0);
		rng = new Random();
		player = new PlayerSkeleton();
		population = new ArrayList<ArrayList<FeatureWeightPair>>(populationSize);
		fitnessResults = new ArrayList<FitnessAssessment>(populationSize);

		for (int i = 0; i < populationSize; i++) {
			ArrayList<FeatureWeightPair> individual = generateRandomIndividual();
			// normalize(individual);
			population.add(individual);
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

		// individual.add(new FeatureWeightPair(new PlayerSkeleton.AggHeight(), randomFloat(-5.0f, 0.0f), false));
		// individual.add(new FeatureWeightPair(new PlayerSkeleton.Bumpiness(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.BumpinessSquared(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.MaxHeight(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.MeanHeight(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.NumHoles(), randomFloat(-5.0f, 0.0f), false));
individual.add(new FeatureWeightPair(new PlayerSkeleton.RowsCleared(), randomFloat(0.0f, 5.0f), true)); // this increases happiness
		individual.add(new FeatureWeightPair(new PlayerSkeleton.RowTransitions(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.StdDevHeight(), randomFloat(-5.0f, 0.0f), false));
individual.add(new FeatureWeightPair(new PlayerSkeleton.SumOfPitDepth(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.TotalHoleDepth(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.WellSums(), randomFloat(-5.0f, 0.0f), false)); // i think this is a superset of sum of pit depth

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

	// returns a new random float according to the gaussian distribution with the configured mu and sigma
	private float nextGaussian(float mu, float sigma) {
		return mu + ((float) rng.nextGaussian() * sigma);
	}

	// run the genetic algorithm for a specified number of generations
	// it will terminate early if convergence is detected. 
	// ConvergenceThreshhold controls how many generations of no better solution being found before terminating
	// returns fitness information for the best individual in the last generation
	public FitnessAssessment trainFor(int generations, int convergenceThreshhold) {
		assert(generations > 0 && convergenceThreshhold > 0);
		System.out.println("training for " + generations + " generations with convergence threshhold of " + convergenceThreshhold + " generations:");
		log("training for " + generations + " generations with convergence threshhold of " + convergenceThreshhold + " generations:");
		// keep track of convergence to an optimal solution
		ArrayList<FeatureWeightPair> bestSoFar = null;
		int generationBestWasFound = 0;
		boolean convergence = false;

		for (int generation = 0; generation < generations && convergence == false; generation++) {
			System.out.println("currently on generation " + generation);
			log("currently on generation " + generation);

			// compute the fitness of everyone
			assessFitnessOfPopulation();
			/*
			fitnessResults.clear(); // clear the results from previous population

			for (ArrayList<FeatureWeightPair> individual : population) {
				rng.setSeed(seed);
				fitnessResults.add(assessFitness(individual));
			}
			 */		
			Collections.sort(fitnessResults);

			log("fitness scores for this generation:");
			for (int j=0; j<fitnessResults.size(); j++)
				log(fitnessResults.get(j).toString());

			System.out.println("The best individual this generation is ");
			System.out.println(fitnessResults.get(fitnessResults.size() - 1));

			ArrayList<FeatureWeightPair> curBest = fitnessResults.get(fitnessResults.size() - 1).individual; // best for this generation
			if (curBest != bestSoFar) { // we just need equality by address
				bestSoFar = curBest;
				generationBestWasFound = generation;
			}

			convergence = generation - generationBestWasFound >= convergenceThreshhold;	

			if (convergence == false) {
				population = reproduce();
				log("population after reproduction:");
				logPopulation();
			}
		}

		if (convergence) {
			System.out.println("Training stopped because convergence was detected with no improvement after " + convergenceThreshhold + " generations");
			log("Training stopped because convergence was detected with no improvement after " + convergenceThreshhold + " generations");
		}
		service.shutdown();
		return fitnessResults.get(fitnessResults.size() - 1);
	}

	// reproduces children
	private ArrayList<ArrayList<FeatureWeightPair>> reproduce() {
		// after the fitnessResults array is sorted, copy the best numElites individuals
		int iterations = (population.size() - numElites) / 2; // how many we need to produce replacements

		ArrayList<ArrayList<FeatureWeightPair>> children = new ArrayList<ArrayList<FeatureWeightPair>>(population.size()); // the next generation

		// copy numElites which should be at the end of the array after sorting
		for (int i = 0; i < numElites; i++)
			children.add(fitnessResults.get(fitnessResults.size() - 1 - i).individual);

		for (int i = 0; i<iterations; i++) {
			// find 2 parents to mate
			ArrayList<FeatureWeightPair> child1 = deepCopyIndividual(tournamentSelection(tournamentSize).individual);
			ArrayList<FeatureWeightPair> child2 = deepCopyIndividual(tournamentSelection(tournamentSize).individual);

			if (isLogging)
				log("mating " + getIndividualAsStr(child1) + " with " + getIndividualAsStr(child2));

			uniformCrossover(child1, child2);

			if (isLogging)
				log("After crossover: child1 = " + getIndividualAsStr(child1) + ", child2 = " + getIndividualAsStr(child2));

			mutate(child1);
			mutate(child2);
			if (isLogging)
				log("After mutation: child1 = " + getIndividualAsStr(child1) + ", child2 = " + getIndividualAsStr(child2));

// normalize(child1);
// normalize(child2);
// if (isLogging)
// log("After normalization: child1 = " + getIndividualAsStr(child1) + ", child2 = " + getIndividualAsStr(child2));
			
			children.add(child1);
			children.add(child2);
		}

		return children;
	}

	// does a deep copy of an individual's weights
	ArrayList<FeatureWeightPair> deepCopyIndividual(ArrayList<FeatureWeightPair> individual) {
		ArrayList<FeatureWeightPair> copy = new ArrayList<FeatureWeightPair>(individual.size());

		for (int i = 0; i < individual.size(); i++)
			copy.add(new FeatureWeightPair(individual.get(i).feature, individual.get(i).weight, individual.get(i).increasesHappiness));

		return copy;
	}

	// returns information about the lowest, average and highest score on an individual after playing numGames games, each game with random piece sequences
	private FitnessAssessment assessFitness(ArrayList<FeatureWeightPair> individual) {
		ArrayList<Integer> scores = new ArrayList<>(numGames);
		ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(numGames);

		for (int game=0; game<numGames; game++) { 
			PlayerSkeleton player = new PlayerSkeleton();
			player.setFeatureWeightPairs(individual);
			player.setSeed(rng.nextLong()); // this is the actual seed used in the sequence
			tasks.add(service.submit(player));
		}

		for (int game=0; game<numGames; game++)
			try {
				scores.add(tasks.get(game).get());
			} catch (Exception e) {
				log("the future computing a game's score was interupted");
				System.out.println("the future computing a game's score was interupted");
			}

		int lowest = scores.get(0);
		int highest = lowest;
		int sum = 0;

		for (Integer score: scores) {
			if (score > highest) {
				highest = score;
			}
			if (score < lowest) {
				lowest = score;
			}
			sum += score;
		}

		return new FitnessAssessment(individual, lowest, sum/numGames, highest);	
	}

	// returns information about the lowest, average and highest score on the entire population after playing numGames games, each game with random piece sequences
	private void assessFitnessOfPopulation() {
		fitnessResults.clear();

		ArrayList<Integer> scores = new ArrayList<>(numGames * population.size());
		ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(numGames);
		long[] seeds = new long[numGames];

		for (int game=0; game<numGames; game++)
			seeds[game] = rng.nextLong();

		for (ArrayList<FeatureWeightPair> individual : population) {
			for (int game=0; game<numGames; game++) { 
				PlayerSkeleton player = new PlayerSkeleton();
				player.setFeatureWeightPairs(individual);
				player.setSeed(seeds[game]); 
				tasks.add(service.submit(player));
			}
		}

		for (int i=0; i<population.size(); i++) {
			for (int game=0; game<numGames; game++)
				try {
					int index = (i * numGames) + game;
					scores.add(tasks.get(index).get());
				} catch (Exception e) {
					log("the future computing a game's score was interupted");
					System.out.println("the future computing a game's score was interupted");
				}
		}

		for (int individual=0; individual<population.size(); individual++) {
			int startIndex = individual * numGames; // the individual's scores are in indices [startIndex, startIndex+numGames)
			int lowest = scores.get(startIndex); // the first game played by this individual
			int highest = lowest;
			int sum = lowest;
			// check second game and up
			for (int game=1; game< numGames; game++) {
				int score = scores.get(startIndex+game);
				if (score > highest) {
					highest = score;
				}
				if (score < lowest) {
					lowest = score;
				}
				sum += score;
			}

			fitnessResults.add(new FitnessAssessment(population.get(individual), lowest, sum/numGames, highest));	
		}

	}

	// returns the fitness assessment of the individual being selected through tournament selection
	private FitnessAssessment tournamentSelection(int tournamentSize) {
		int best = randomInt(0, population.size());
		FitnessAssessment bestFitness = fitnessResults.get(best);

		for (int i = 2; i <= tournamentSize; i++) {
			int next = randomInt(fitnessResults.size()-numElites, fitnessResults.size());

			// ensure that the next individual selected is different from the existing best
			while (best == next)
				next = randomInt(fitnessResults.size()-numElites, fitnessResults.size());

			FitnessAssessment fitness = fitnessResults.get(next);

			if (fitness.compareTo(bestFitness)  > 0) {
				bestFitness = fitness;
				best = next;
			}
		}

		return bestFitness;
	}

	// crosses over 2 individuals using uniform crossover
	private void uniformCrossover(ArrayList<FeatureWeightPair> x, ArrayList<FeatureWeightPair> y) {
		for (int i = 0; i < x.size(); i++) {
			if (crossoverRate >= randomFloat(0.0f, 1.0f)) {
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
			if (individual.get(i).increasesHappiness)
				do {
					n  = nextGaussian(0.0f, mutationSigma);
				} while (individual.get(i).weight + n <= 0.0f);

			else
				do {
					n  = nextGaussian(0.0f, mutationSigma);
				} while (individual.get(i).weight + n >= 0.0f);

			individual.get(i).weight += n;
		}
	}

	// computes the length of a weight vector
	private static double length(ArrayList<FeatureWeightPair> vec) {
		float sumSquares = 0.0f;
		for (FeatureWeightPair f : vec) 
			sumSquares += f.weight * f.weight;

		return Math.sqrt(sumSquares);
	}

	// normalizes a vector
	public static void normalize(ArrayList<FeatureWeightPair> vec) {
		double length = length(vec);
		
		for (FeatureWeightPair f : vec)
			f.weight /= length;
	}

	public static void loggerInit(String filename) {
		try {
			FileHandler handler = new FileHandler(filename);
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
		int cores = Runtime.getRuntime().availableProcessors();
		System.out.println("Running genetic algorithm on a machine with " + cores + " logical cores:");
		Scanner sc = new Scanner(System.in);
		String filename;
		System.out.print("\nEnter log filename: ");
		filename = sc.nextLine();
		loggerInit(filename);
		int elites = 0, games = 0, populationSize = 0, tournamentSize = 0, generations = 0, convergenceThreshhold = 0;
		float mutationSigma = 0.0f, crossoverRate = 0.0f;
		System.out.println("Enter parameters");
		System.out.print("\nnumber of generations: ");
		generations = sc.nextInt();
		System.out.print("\nnumber of generations with no improvement  before training stops (convergence): ");
		convergenceThreshhold = sc.nextInt();
		System.out.print("\nCrossover rate: ");
		crossoverRate = sc.nextFloat();
		System.out.print("\nNumber of elites (population size - elites should be even): ");
		elites = sc.nextInt();
		System.out.print("\nnumber of games per individual: ");
		games = sc.nextInt();
		System.out.print("\nstandard deviation for mutation: ");
		mutationSigma = sc.nextFloat();
		System.out.print("Population size: ");
		populationSize = sc.nextInt();
		System.out.print("\nTournament size: ");
		tournamentSize = sc.nextInt();

		GeneticAlgorithm ga = new GeneticAlgorithm(crossoverRate, elites, games, mutationSigma, populationSize, tournamentSize);
		long startTime = System.currentTimeMillis();
		FitnessAssessment result =ga.trainFor(generations, convergenceThreshhold);
		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("Training complete. The best individual is ");
		System.out.println(result);
		System.out.println("GA took: " + totalTime + " ms (" + totalTime/1000 + " seconds)");
		sc.close();
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
			if (average < other.average)
				return -1;

			if (average > other.average)
				return 1;

			// tiebreak using lowest
			if (lowest < other.lowest)
				return -1;

			if (lowest > other.lowest)
				return 1;

			// then highest
			if (highest < other.highest)
				return -1;

			if (highest > other.highest)
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
