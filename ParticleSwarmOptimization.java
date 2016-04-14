// main reference: essentials of metaheuristics, page 59
// this implementation is essentially straight from the book minus informants and particle jump size
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

public class ParticleSwarmOptimization {
	public static final boolean isLogging = true;
	private static final Logger logger = Logger.getLogger("PSO");
	// since particles are just vectors in r^n, the correspondence between the i^{th} coordinate 
	// and the corresponding features and its min and max allowed weights is maintained here
	private ArrayList<FeatureBoundaryPair> features;
	private ArrayList<Particle> population;
	private int populationSize; // size of the swarm
	// inertia is the proportion of velocity to be retained; or its "momentum" so that the particle's movement isn't too erratic 
	// Low settings favour exploitation (local search) while a high value favours global search.
	// A decaying inertia weight can be used, with the aim of favoring global search at the start of the algorithm and local search later.
	// If inertia is not reduced with time, choose a value in [0.8, 1.2]
	// see http://tracer.uc3m.es/tws/pso/parameters.html
	private float startInertia, endInertia, inertia;
	// the cognitive component is the proportion of the particle's personal best to be retained 
	// allows a particle's memory of where it has previously found good solutions to influence its movement direction
	// If this is large, particles tend to move more towards their own personal bests rather than towards global bests.
	// This breaks the swarm into a lot of separate hill-climbers rather than a joint searcher.
	private float cognitive;
	// the social component is the proportion of global best to be retained.
	// allows the knowledge of other swarm members (i.e. where other members of the swarm have found good solutions) to influence a particle's movement
	// If this is large, particles tend to move more towards the best known region. This converts the algorithm into one large hill-climber
	// the book suggests setting it to 0, other sources disagree on this; we'll have to experiment
	private float social;
	// other sources on parameter selection:
	// http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=4019136&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D4019136
// 		http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=5957817&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D5957817
		// http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=870279&tag=1&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D870279%26tag%3D1

	public int numGames = 0; // number of games to run to assess fitness of an individual
	private Random rng;
	private ArrayList<FitnessAssessment> fitnessResults;
	private ExecutorService service = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
	
	ParticleSwarmOptimization(int populationSize, int games, float cognitive, float startInertia, float endInertia, float social) {
		rng = new Random();
		this.populationSize = populationSize;
		this.numGames = games;
		this.cognitive = cognitive;
		this.startInertia= startInertia;
		this.endInertia= endInertia;
		inertia = startInertia;
		this.social = social;
		population = new ArrayList<Particle>(populationSize);
		fitnessResults = new ArrayList<FitnessAssessment>(populationSize);
		features = new ArrayList<FeatureBoundaryPair>();
		// add features here
		// features.add(new FeatureBoundaryPair(new PlayerSkeleton.Bumpiness(), -10.0f, 0.0f));
		features.add(new FeatureBoundaryPair(new PlayerSkeleton.BumpinessSquared(), -10.0f, 0.0f));
		features.add(new FeatureBoundaryPair(new PlayerSkeleton.MaxHeight(), -10.0f, 0.0f));
		features.add(new FeatureBoundaryPair(new PlayerSkeleton.MeanHeight(), -10.0f, 0.0f));
		features.add(new FeatureBoundaryPair(new PlayerSkeleton.NumHoles(), -10.0f, 0.0f));
features.add(new FeatureBoundaryPair(new PlayerSkeleton.RowsCleared(), 0.0f, 10.0f));
		features.add(new FeatureBoundaryPair(new PlayerSkeleton.RowTransitions(), -10.0f, 0.0f));
		features.add(new FeatureBoundaryPair(new PlayerSkeleton.StdDevHeight(), -10.0f, 0.0f));
		features.add(new FeatureBoundaryPair(new PlayerSkeleton.SumOfPitDepth(), -10.0f, 0.0f));
		features.add(new FeatureBoundaryPair(new PlayerSkeleton.TotalHoleDepth(), -10.0f, 0.0f));
features.add(new FeatureBoundaryPair(new PlayerSkeleton.WellSums(), -10.0f, 0.0f));

		generateRandomPopulation();
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
	
	private String getIndividualAsStr(Particle individual) {
		String str = "";
		
		for (int i=0; i<individual.position.length; i++)
			str += individual.position[i] + ", ";

		str += " velocity: ";
		
		for (int i=0; i<individual.position.length; i++)
			str += individual.velocity[i] + ", ";

str += " personal best: ";

for (int i=0; i<individual.position.length; i++)
	str += individual.personalBestVector[i] + ", ";

str += " personal best score: " + individual.personalBestFitness;

		return str;
	}
	
	private String getIndividualAsStr(int i) {
		return getIndividualAsStr(population.get(i));
	}
	
	// returns a new float in the range [min, max)
		private float randomFloat(float minInclusive, float maxExclusive) {
			return (float) ThreadLocalRandom.current().nextDouble(minInclusive, maxExclusive);
		}
		
		// generate the initial particle swarm
		private void generateRandomPopulation() {
			for (int i=0; i<populationSize; i++)
				population.add(generateRandomStationaryParticle());
		}
		
		// generates a random particle with a 0 velocity vector
		private Particle generateRandomStationaryParticle() {
			float[] pos = new float[features.size()];
			
			for (int i=0; i<features.size(); i++)
				pos[i] = randomFloat(features.get(i).min, features.get(i).max);
			
			// normalize(pos);
			return new Particle(pos);
		}
		
		// computes the length of a weight vector
		private static double length(float[] vec) {
			float sumSquares = 0.0f;
			for (float f : vec) 
				sumSquares += f;

			return Math.sqrt(sumSquares);
		}

		// normalizes a vector
		public static void normalize(float[] vec) {
			double length = length(vec);
			
			for (int i=0; i<vec.length; i++)
				vec[i] /=  length;
		}


		// run the algorithm for a specified number of generations
public FitnessAssessment trainFor(int generations, int convergenceThreshhold) {
	System.out.println("Training for " + generations + " generations with convergence threshhold of " + convergenceThreshhold + " generations:");
	log("training for " + generations + " generations with convergence threshhold of " + convergenceThreshhold + " generations:");
	
	float[] bestPosition = population.get(0).getCopyOfPosition();
	// initialize to the worst possible score
	FitnessAssessment bestFitness = new FitnessAssessment(bestPosition, 0, 0, 0);
	int generationBestWasFound = 0;
	boolean convergence = false;
	float changeInInertia = (endInertia - startInertia) / generations; 
	
	for (int generation = 0; generation < generations && convergence == false; generation++) {
		System.out.println("currently on generation " + generation);
		long[] seeds = new long[numGames];
		
		for (int game=0; game<numGames; game++)
			seeds[game] = rng.nextLong();
		
		 bestFitness = assessFitness(seeds, bestPosition);
		assessFitnessOfPopulation(seeds);
		
		for (int i=0; i<population.size(); i++) {
			Particle particle = population.get(i);
			FitnessAssessment curFitness = fitnessResults.get(i);
			if (curFitness.compareTo(bestFitness) > 0) {
				bestFitness = curFitness;
				bestPosition = particle.getCopyOfPosition();
				generationBestWasFound = generation;
			}
			
			if (curFitness.compareTo(particle.personalBestFitness) > 0) {
				particle.personalBestFitness = curFitness;
				particle.setBestFitnessVector (particle.position); // we need a deep copy so position changes don't affect this
		}
	}
		
		System.out.println("The best particle found so far is: ");
		System.out.println(bestFitness);

		convergence = generation - generationBestWasFound >= convergenceThreshhold;
	
		if (convergence == false) {
		determineVelocity(bestPosition);
		mutate();
		log("population after mutation:");
		logPopulation();
		}
		inertia += changeInInertia;

	}
	
	if (convergence) {
		System.out.println("Training stopped because convergence was detected with no improvement after " + convergenceThreshhold + " generations");
	log("Training stopped because convergence was detected with no improvement after " + convergenceThreshhold + " generations");
	}
	service.shutdown();
	
	return bestFitness;
}
		
// assess fitness of an individual
private FitnessAssessment assessFitness(long[] seeds, float[] pos) {
	ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(numGames);
				ArrayList<Integer> scores = new ArrayList<Integer>(numGames);

				for (int game=0; game<numGames; game++) { 
					PlayerSkeleton player = getNewPlayerUsingWeights(pos);
					player.setSeed(seeds[game]); 
					tasks.add(service.submit(player));
				}

	// retrieve scores for best vector
				for (int game=0; game<numGames; game++)
					try {
					scores.add(tasks.get(game).get());
							} catch (Exception e) {
								log("the future computing a game's score was interupted");
								System.out.println("the future computing a game's score was interupted");
				}
				
	int lowest = scores.get(0); // the first game played by best vector
				int highest = lowest;
				int sum = lowest;
	//check second game and up
							for (int game=1; game< numGames; game++) {
		int score = scores.get(game);
					if (score > highest) {
						highest = score;
					}
					if (score < lowest) {
						lowest = score;
					}
					sum += score;
				}

	return new FitnessAssessment(pos, lowest, sum/numGames, highest);	
}


// returns a PlayerSkeleton initialized with a weight vector
private PlayerSkeleton getNewPlayerUsingWeights(float[] weights) {
	PlayerSkeleton player = new PlayerSkeleton();
	int numFeatures = features.size();
	ArrayList<FeatureWeightPair> fwPairs = new ArrayList<FeatureWeightPair>(numFeatures);
	
	for (int i=0; i<numFeatures; i++)
	fwPairs.add(new FeatureWeightPair(features.get(i).feature, weights[i], false)); // the third argument doesn't really matter here; only used by GA
	
		player.setFeatureWeightPairs(fwPairs);
	return player;
}

// computes the fitness of all particles
		private void assessFitnessOfPopulation(long[] seeds) {
			fitnessResults.clear();
		
			ArrayList<Integer> scores = new ArrayList<>(numGames * population.size());
			ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(numGames * population.size());

			for (Particle individual : population) {
			for (int game=0; game<numGames; game++) { 
				PlayerSkeleton player = getNewPlayerUsingWeights(individual.position);
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
//check second game and up
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

			fitnessResults.add(new FitnessAssessment(population.get(individual).position, lowest, sum/numGames, highest));	
		}
			
			tasks.clear();
			scores.clear();
			// update the score of the personal best positions of the particle
			// we can't reuse the old value as the personal best would perform differently across generations
			for (Particle individual : population) {
			for (int game=0; game<numGames; game++) { 
				PlayerSkeleton player = getNewPlayerUsingWeights(individual.personalBestVector);
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
//check second game and up
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
			
			FitnessAssessment f = new FitnessAssessment(population.get(individual).personalBestVector, lowest, sum/numGames, highest);
			population.get(individual).personalBestFitness = f;
		}
			
		}

		// determines the velocity of particles in this generation
		// bestSoFar is the coordinates of best particle found so far
		private void determineVelocity(float[] bestSoFar) {
			for (Particle p : population)
				for (int i=0; i<features.size(); i++) {
					float c1 = randomFloat(0.0f, cognitive);
					float c2 = randomFloat(0.0f, social);
					float diffFromPersonalBest = p.personalBestVector[i] - p.position[i];
					float diffFromBestFoundSoFar = bestSoFar[i] - p.position[i];
					p.velocity[i] = (inertia * p.velocity[i]) + (c1 * diffFromPersonalBest) + (c2 * diffFromBestFoundSoFar);
					
				}
		}
		
		// mutates the particles by applying their velocities to their current position
		private void mutate() {
			for (Particle p : population)
				for (int i=0; i<features.size(); i++) {
					float newCoord = p.position[i] + p.velocity[i];
					newCoord = clamp(newCoord, features.get(i).min, features.get(i).max); // prevent particles from escaping boundaries
					p.position[i] = newCoord;
		}
		}
		
		// clamps a value if it is outside [min, max] so that it falls inside this range
		private static float clamp(float val, float min, float max) {
		    return Math.max(min, Math.min(max, val));
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
			System.out.println("Running particle swarm optimization on a machine with " + cores + " logical cores:");
			Scanner sc = new Scanner(System.in);
			
			int games = 0, populationSize = 0, generations = 0, convergenceThreshhold = 0;
			String filename;
			float cognitive = 0.0f, startInertia = 0.0f, endInertia = 0.0f, social = 0.0f;
	        System.out.println("Enter parameters");
	        System.out.print("\nlog file name: ");
	        filename = sc.nextLine();
	        System.out.print("\nnumber of generations: ");
	        generations = sc.nextInt();
	        System.out.print("\nnumber of generations with no improvement  before training stops (convergence): ");
	        convergenceThreshhold = sc.nextInt();
	        System.out.print("\nnumber of games per individual: ");
	        games = sc.nextInt();
	        System.out.print("Population size: ");
	        populationSize = sc.nextInt();
	        System.out.print("\ncognitive (>= 0.0): ");
	        cognitive = sc.nextFloat();
	        System.out.print("\nstarting inertia (>= 0.0): ");
	        startInertia = sc.nextFloat();
	        System.out.print("\nending inertia (>= 0.0): ");
	        endInertia = sc.nextFloat();
			
			System.out.print("\nsocial (>= 0.0): ");
			social = sc.nextFloat();
			
			loggerInit(filename);
			ParticleSwarmOptimization pso = new ParticleSwarmOptimization(populationSize, games, cognitive, startInertia, endInertia, social);
			long startTime = System.currentTimeMillis();
			FitnessAssessment result =pso.trainFor(generations, convergenceThreshhold);
	        long endTime   = System.currentTimeMillis();
	        long totalTime = endTime - startTime;
			System.out.println("Training complete. The best individual is ");
			System.out.println(result);
	        System.out.println("PSO took: " + totalTime + " ms (" + totalTime/1000 + " seconds)");
	        sc.close();
		}

	// a utility class that associates a feature with the minimum and maximum weight
	// this constrains the search space
	private class FeatureBoundaryPair {
		float min = 0.0f, max = 0.0f;
		FeatureFunction feature;
		
		FeatureBoundaryPair(FeatureFunction f, float min, float max) {
			assert(min <= max);
			this.feature = f;
			this.min = min;
			this.max = max;
		}
	}

	// stores information about the fitness of an individual
		public static class FitnessAssessment implements Comparable<FitnessAssessment> {
			public float[] individual;
			public int lowest, highest; // scores
			public float average; // average score

			public FitnessAssessment(float[] individual, int l, float a, int h) {
				this.individual = individual;
				lowest = l;
				average = a;
				highest = h;
			}

			public int compareTo(ParticleSwarmOptimization.FitnessAssessment other) {
				
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

				for (float w : individual)
					str += w + ", ";

				return str;
			}
		}

}
