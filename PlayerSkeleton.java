import java.util.*;
import java.util.concurrent.Callable;

public class PlayerSkeleton implements Callable<Integer> {

	public static final int COLS = 10;
	public static final int ROWS = 21;
	public static final int N_PIECES = 7;
	private ArrayList<FeatureWeightPair> features;
	private Long seed;

	public PlayerSkeleton() {
		features = new ArrayList<FeatureWeightPair>();
	}

	public void setFeatureWeightPairs(ArrayList<FeatureWeightPair> features) {
		this.features = features;
	}

	// returns h(n) for the given state
	private float evaluate(ImprovedState s) {
		if (s.hasLost())
			return Float.NEGATIVE_INFINITY;

		float sum = 0.0f;

		for (FeatureWeightPair f : features) {
			sum += f.feature.evaluate(s) * f.weight;
		}
		return sum;
	}

	public Integer call() throws InterruptedException {
		return playGameWithImprovedState();
	}

	// lookahead is the number of future pieces we will consider
	// for example, a value of 1 means that for the current state, legal move m
	// is executed to get some new state.
	// From this new state, calculate the highest utility each of the 7 pieces
	// can have and average this out to get the utility for move m
	private MoveUtilityPair pickMove(State s, int nextPiece, int lookahead) {
		return pickMove(new ImprovedState(s), nextPiece, lookahead);
	}

	private MoveUtilityPair pickMove(ImprovedState s, int nextPiece, int lookahead) {
		float best = Float.NEGATIVE_INFINITY;
		int[][][] movesForAllPieces = s.getLegalMovesForAllPieces();
		int bestMove = 0;

		if (lookahead == 0) {
			for (int move = 0; move < movesForAllPieces[nextPiece].length; move++) {
				s.makeMove(movesForAllPieces[nextPiece][move]);
				float utility = evaluate(s);
				s.undo();

				if (utility > best) {
					best = utility;
					bestMove = move;
				}
			}
		} else {
			// evaluate the score of the state resulting from each possible move

			ArrayList<MoveUtilityPair> scores = new ArrayList<MoveUtilityPair>(movesForAllPieces[nextPiece].length);
			int numNodesToExpand = 3; // the best n nodes to expand

			for (int move = 0; move < movesForAllPieces[nextPiece].length; move++) {
				s.makeMove(movesForAllPieces[nextPiece][move]);
				float score = evaluate(s);
				scores.add(new MoveUtilityPair(move, score));
				s.undo();
			}

			Collections.sort(scores);
			// indices lower than this not expanded
			int minIndex = (scores.size() - numNodesToExpand >= 0) ? scores.size() - numNodesToExpand : 0;

			for (int i = minIndex; i < scores.size(); i++) {
				int move = scores.get(i).move;
				float curUtility = scores.get(i).utility;
				if (curUtility > Float.NEGATIVE_INFINITY) { // no point if we
															// have already lost

					ImprovedState future = s.tryMove(movesForAllPieces[nextPiece][move]);
					float sum = 0.0f; // sum of future utilities

					for (int piece = 0; piece < ImprovedState.N_PIECES; piece++) {
						future.setNextPiece(piece);
						sum += pickMove(future, piece, lookahead - 1).utility;
					}
					float avgFutureUtility = sum / ImprovedState.N_PIECES;

					if (avgFutureUtility > best) {
						best = avgFutureUtility;
						bestMove = move;
					}
				}
			}
		}

		return new MoveUtilityPair(bestMove, best);
	}

	// plays a game , returning the number of rows complete
	// use the setFeatureWeightPairs function first
	public int playGame(boolean alwaysDraw, boolean drawOnLoss) throws InterruptedException {
		assert (!features.isEmpty()); // must set some features to use first
		State s = new State();
		if (alwaysDraw)
			new TFrame(s);

		while (!s.hasLost()) {
			s.makeMove(pickMove(s, s.getNextPiece(), 0).move);
			if (alwaysDraw) {
				s.draw();
				s.drawNext(0, 0);

				try {
					Thread.sleep(100);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		if (alwaysDraw == false && drawOnLoss == true) {
			new TFrame(s);
			s.draw();
			s.drawNext(0, 0);
		}

		return s.getRowsCleared();
	}

	// sets the seed for playing a game with ImprovedState
	public void setSeed(long seed) {
		this.seed = seed;
	}

	// identical to the playGame function except that ImprovedState is used
	// which is much more useful for training purposes
	// note that drawing isn't supported though with this one
	// use the setFeatureWeightPairs function first
	public int playGameWithImprovedState() throws InterruptedException {
		assert (!features.isEmpty());
		ImprovedState s = new ImprovedState();
		if (seed != null) {
			s.setSeed(seed);
		}
		// s.setSeed(1459523385737L);
		s.pickNextPiece();

		while (!s.hasLost()) {
			s.makeMove(pickMove(s, s.getNextPiece(), 0).move);
			s.pickNextPiece();
		}

		return s.getRowsCleared();
	}

	public static void main(String[] args) throws InterruptedException {
		PlayerSkeleton p = new PlayerSkeleton();

		ArrayList<FeatureWeightPair> fwPairs = new ArrayList<FeatureWeightPair>();

		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.MaxHeight(), -0.87448764f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.MeanHeight(), -0.5054159f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.NumHoles(), -7.5740294f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.RowTransitions(), -1.1404463f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.StdDevHeight(), -1.9600217f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.SumOfPitDepth(), -1.4352853f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.WellSums(), -0.9274553f, false));
		p.setFeatureWeightPairs(fwPairs);

		long startTime = System.currentTimeMillis();
		System.out.println("You have completed " + p.playGameWithImprovedState() + " rows.");
		// System.out.println("You have completed "+p.playGame(false, true) +"
		// rows.");
		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("PlayerSkeleton took: " + totalTime + "ms");

	}

	// Considered as a pit if the adjacent columns are >= 2. Depth = diff in
	// height with the shortest adjacent col
	// should refactor this some more
	public static class SumOfPitDepth implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int[] top = s.getTop();
			int sumOfPitDepth = 0;

			int pitHeight;
			int heightOfLeftColumn;
			int heightOfRightCol;

			// pit depth of first column
			pitHeight = top[0];
			heightOfRightCol = top[1];
			int heightDiff = heightOfRightCol - pitHeight;
			if (heightDiff > 2) {
				sumOfPitDepth += heightDiff;
			}

			for (int col = 0; col < State.COLS - 2; col++) {
				heightOfLeftColumn = top[col];
				pitHeight = top[col + 1];
				heightOfRightCol = top[col + 2];

				int leftDiff = heightOfLeftColumn - pitHeight;
				int rightDiff = heightOfRightCol - pitHeight;
				int minDiff = Math.min(leftDiff, rightDiff);

				if (minDiff >= 2) {
					sumOfPitDepth += minDiff;
				}
			}

			// pit depth of last column
			pitHeight = top[State.COLS - 1];
			heightOfLeftColumn = top[State.COLS - 2];
			heightDiff = heightOfLeftColumn - pitHeight;
			if (heightDiff > 2) {
				sumOfPitDepth += heightDiff;
			}

			return (float) sumOfPitDepth;
		}
	}

	// Returns the average height across the columns
	public static class MeanHeight implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			float avgHeight = 0;
			int[] top = s.getTop();
			for (int i = 0; i < State.COLS; i++) {
				avgHeight += top[i];
			}
			return avgHeight / State.COLS;
		}
	}

	public static class RowsCleared implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			return (float) s.getRowsCleared();
		}
	}

	// returns aggregate height for all columns
	public static class AggHeight implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int aggHeight = 0;

			int[] top = s.getTop();

			for (int i = 0; i < State.COLS; i++) {
				aggHeight += top[i];
			}

			return (float) aggHeight;
		}
	}

	// returns number of holes. A hole is an empty space such that there is at
	// least one tile in the same column above it
	public static class NumHoles implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int numHoles = 0;

			int[][] field = s.getField();
			int[] top = s.getTop();

			for (int c = 0; c < State.COLS; c++) {
				for (int r = 0; r < top[c]; r++) {
					if (field[r][c] == 0) {
						numHoles++;
					}
				}
			}

			return numHoles;
		}
	}

	// number of occupied cells above holes
	public static class TotalHoleDepth implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int holeDepths = 0;

			int[][] field = s.getField();
			int[] top = s.getTop();

			for (int c = 0; c < State.COLS; c++) {
				for (int r = 0; r < top[c]; r++)
					if (field[r][c] == 0) { // this is a hole, check rows above
						for (int i = r + 1; i <= top[c]; i++)
							if (field[i][c] != 0)
								holeDepths++;
					}
			}

			return holeDepths;
		}
	}

	// calculates bumpiness, the sum of the absolute differences between heights
	// of consecutive adjacent columns
	public static class Bumpiness implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int bumpiness = 0;

			int[] top = s.getTop();

			for (int i = 0; i < State.COLS - 1; i++) {
				bumpiness += Math.abs(top[i] - top[i + 1]);
			}

			return bumpiness;
		}
	}

	// The name is a bit misleading, but i can't think of a better name
	// computes sum of the absolute differences between squared heights of
	// consecutive adjacent columns
	// for example, for columns with heights of 5 and 4, the difference would be
	// 5^2 - 4^2 = 9
	// the hope is that this encourages "smoother" play
	public static class BumpinessSquared implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int bumpiness = 0;

			int[] top = s.getTop();

			for (int i = 0; i < State.COLS - 1; i++) {
				bumpiness += Math.abs((top[i] * top[i]) - (top[i + 1] * top[i + 1]));
			}

			return bumpiness;
		}
	}

	// maximum column height
	public static class MaxHeight implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int maxColumnHeight = 0;

			int[] top = s.getTop();

			for (int i = 0; i < State.COLS; i++) {
				if (top[i] > maxColumnHeight) {
					maxColumnHeight = top[i];
				}
			}
			return maxColumnHeight;
		}
	}

	// computes total row transitions. Row transitions happen when an empty cell
	// is adjacent to a filled cell and vice versa
	public static class RowTransitions implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int nRowTransitions = 0;
			int[][] field = s.getField();

			for (int r = 0; r < State.ROWS; r++)
				for (int c = 0; c < State.COLS - 1; c++) {
					boolean isCurEmpty = field[r][c] == 0, isNextEmpty = field[r][c + 1] == 0;

					if ((isCurEmpty && !isNextEmpty) || (!isCurEmpty && isNextEmpty))
						nRowTransitions++;
				}
			return nRowTransitions;
		}
	}

	// computes the standard deviation of column heights
	public static class StdDevHeight implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			float sum = 0.0f;
			int[] top = s.getTop();

			for (int i = 0; i < State.COLS; i++)
				sum += top[i];

			float mean = sum / State.COLS;
			float sumDiffFromMean = 0.0f;

			for (int c = 0; c < State.COLS; c++)
				sumDiffFromMean = (top[c] - mean) * (top[c] - mean);

			return (float) Math.sqrt(sumDiffFromMean / (State.COLS - 1));
		}
	}

	// computes total well sums
	// a well is the number of empty cells above a column's top piece such that
	// the top cell in the sequence is surrounded by either the board boundary
	// or occupied cells in neighbouring columns
	public static class WellSums implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int wellSums = 0;
			int[] top = s.getTop();
			// check well sums for columns that don't boarder the left and right
			// edge
			for (int c = 1; c < State.COLS - 1; c++) {
				int curTop = top[c];
				// get the min height of the adjacent left and right columns
				int minAdjColHeight = Math.min(top[c - 1], top[c + 1]);
				// if the current column is <= than the heights of its
				// neighbours, then there is a well
				if (curTop <= minAdjColHeight)
					wellSums += minAdjColHeight - curTop; // this works
															// correctly too if
															// the heights of
															// this columns and
															// its neighbours
															// are the same
			}

			// check for the left and rightmost column
			if (top[1] > top[0])
				wellSums += top[1] - top[0];

			if (top[State.COLS - 2] > top[State.COLS - 1])
				wellSums += top[State.COLS - 2] - top[State.COLS - 1];

			return wellSums;
		}
	}

	// convenience class for associating a move index with its utility
	private class MoveUtilityPair implements Comparable<MoveUtilityPair> {
		public float utility = 0.0f;
		public int move = 0;

		MoveUtilityPair(int m, float u) {
			move = m;
			utility = u;
		}

		public int compareTo(MoveUtilityPair other) {
			if (utility < other.utility)
				return -1;

			if (utility > other.utility)
				return 1;

			return 0;
		}
	}
}

/* Benchmarker.java
//runs PlayerSkeleton for n random games to obtain information about agent performance for report
//the output to the console should be redirected to a .csv file, and opened in excel for analysis/drawing

import java.util.ArrayList;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class Benchmarker {
	public static void main(String[] args) throws InterruptedException {
		int games = 0;
		Scanner sc = new Scanner(System.in);
		games = sc.nextInt();
		int cores = Runtime.getRuntime().availableProcessors();
		ExecutorService service = Executors.newFixedThreadPool(cores);
		ArrayList<Integer> scores = new ArrayList<>(games);
		ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(games);

		System.out.println("game number, rows completed");

		for (int i = 0; i < games; i++) {
			PlayerSkeleton p = new PlayerSkeleton();
			ArrayList<FeatureWeightPair> fwPairs = new ArrayList<FeatureWeightPair>();
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.MaxHeight(), -0.87448764f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.MeanHeight(), -0.5054159f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.NumHoles(), -7.5740294f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.RowTransitions(), -1.1404463f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.StdDevHeight(), -1.9600217f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.SumOfPitDepth(), -1.4352853f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.WellSums(), -0.9274553f, false));
			p.setFeatureWeightPairs(fwPairs);
			tasks.add(service.submit(p));
		}

		for (int i = 0; i < games; i++) {
			try {
				int score = tasks.get(i).get();
				// scores.add(score);
				System.out.println(i + ", " + score);
			} catch (Exception e) {
				System.out.println("the future computing a game's score was interupted");
			}
		}

		sc.close();
		service.shutdown();
	}
}
*/

/* GeneticAlgorithm.java
import java.io.FileWriter;
import java.io.IOException;
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

public class GeneticAlgorithm {
	private static final boolean isLogging = true;
	private static final Logger logger = Logger.getLogger("GeneticAlgorithm");
	public int populationSize = 0;
	public float crossoverRate = 0.0f;
	public int numGames = 0; // number of games to run to assess fitness of an
								// individual
	public float mutationSigma = 0.0f; // standard deviation for mutation
	public int tournamentSize = 2; // 2's the most common setting. 1 is random
									// selection, higher values causes higher
									// selection pressure
	public int numElites = 4; // number of numElites to keep
	ArrayList<ArrayList<FeatureWeightPair>> population;
	ArrayList<FitnessAssessment> fitnessResults;
	PlayerSkeleton player;
	private Random rng;
	private ExecutorService service;

	public GeneticAlgorithm(float crossover, int elites, int games, float mutationSigma, int populationSize,
			int tournamentSize, int threads) {
		service = Executors.newFixedThreadPool(threads);
		this.crossoverRate = crossover;
		this.numElites = elites;
		this.numGames = games;
		this.mutationSigma = mutationSigma;
		this.populationSize = populationSize;
		this.tournamentSize = tournamentSize;

		assert (numElites >= 0 && numElites <= populationSize);
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
		for (int i = 0; i < population.size(); i++) {
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

	static FileWriter writer;

	private ArrayList<FeatureWeightPair> generateRandomIndividual() {
		ArrayList<FeatureWeightPair> individual = new ArrayList<FeatureWeightPair>();
		// all the feature functions we're using so far contribute negatively to
		// happiness and so should be minimized,
		// hence their weights should be negative
		// for example, the presence of holes should decrease happiness

		// individual.add(new FeatureWeightPair(new PlayerSkeleton.AggHeight(),
		// randomFloat(-5.0f, 0.0f), false));
		// individual.add(new FeatureWeightPair(new PlayerSkeleton.Bumpiness(),
		// randomFloat(-5.0f, 0.0f), false));
		// individual.add(new FeatureWeightPair(new
		// PlayerSkeleton.BumpinessSquared(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.MaxHeight(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.MeanHeight(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.NumHoles(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.RowsCleared(), randomFloat(0.0f, 5.0f), true)); // this
																												// increases
																												// happiness
		individual.add(new FeatureWeightPair(new PlayerSkeleton.RowTransitions(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.StdDevHeight(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.SumOfPitDepth(), randomFloat(-5.0f, 0.0f), false));
		// individual.add(new FeatureWeightPair(new
		// PlayerSkeleton.TotalHoleDepth(), randomFloat(-5.0f, 0.0f), false));
		individual.add(new FeatureWeightPair(new PlayerSkeleton.WellSums(), randomFloat(-5.0f, 0.0f), false)); // i
																												// think
																												// this
																												// is
																												// a
																												// superset
																												// of
																												// sum
																												// of
																												// pit
																												// depth

		return individual;
	}

	// returns a new float in the range [min, max)
	private float randomFloat(float minInclusive, float maxExclusive) {
		return (float) ThreadLocalRandom.current().nextDouble(minInclusive, maxExclusive);
	}

	// returns a new int in the range [min, max)
	// to have both endpoints included, pass a value of max+1 to this function
	private int randomInt(int minInclusive, int maxExclusive) {
		return ThreadLocalRandom.current().nextInt(minInclusive, maxExclusive);
	}

	// returns a new random float according to the gaussian distribution with
	// the configured mu and sigma
	private float nextGaussian(float mu, float sigma) {
		return mu + ((float) rng.nextGaussian() * sigma);
	}

	// run the genetic algorithm for a specified number of generations
	// it will terminate early if convergence is detected.
	// ConvergenceThreshhold controls how many generations of no better solution
	// being found before terminating
	// returns fitness information for the best individual in the last
	// generation
	public FitnessAssessment trainFor(int generations, int convergenceThreshhold) {
		assert (generations > 0 && convergenceThreshhold > 0);
		System.out.println("training for " + generations + " generations with convergence threshhold of "
				+ convergenceThreshhold + " generations:");
		log("training for " + generations + " generations with convergence threshhold of " + convergenceThreshhold
				+ " generations:");
		// keep track of convergence to an optimal solution
		ArrayList<FeatureWeightPair> bestSoFar = null;
		int generationBestWasFound = 0;
		boolean convergence = false;

		for (int generation = 0; generation < generations && convergence == false; generation++) {
			System.out.println("currently on generation " + generation);
			log("currently on generation " + generation);

			// compute the fitness of everyone
			assessFitnessOfPopulation();
			Collections.sort(fitnessResults);

			log("fitness scores for this generation:");
			for (int j = 0; j < fitnessResults.size(); j++) {
				log(fitnessResults.get(j).toString());
				// if (generation == 0) {
				// writeToCSV(fitnessResults.get(j));
				// }
			}

			System.out.println("The best individual this generation is ");
			System.out.println(fitnessResults.get(fitnessResults.size() - 1));

			ArrayList<FeatureWeightPair> curBest = fitnessResults.get(fitnessResults.size() - 1).individual; // best
																												// for
																												// this
																												// generation
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
			System.out.println("Training stopped because convergence was detected with no improvement after "
					+ convergenceThreshhold + " generations");
			log("Training stopped because convergence was detected with no improvement after " + convergenceThreshhold
					+ " generations");
		}
		service.shutdown();
		return fitnessResults.get(fitnessResults.size() - 1);
	}

	// reproduces children
	private ArrayList<ArrayList<FeatureWeightPair>> reproduce() {

		// after the fitnessResults array is sorted, copy the best numElites
		// individuals
		int iterations = (population.size() - numElites) / 2; // how many we
																// need to
																// produce
																// replacements

		ArrayList<ArrayList<FeatureWeightPair>> children = new ArrayList<ArrayList<FeatureWeightPair>>(
				population.size()); // the next generation

		// copy numElites which should be at the end of the array after sorting
		for (int i = 0; i < numElites; i++) {
			children.add(fitnessResults.get(fitnessResults.size() - 1 - i).individual);
		}

		for (int i = 0; i < iterations; i++) {

			// find 2 parents to mate
			ArrayList<FeatureWeightPair> child1 = deepCopyIndividual(tournamentSelection(tournamentSize).individual);

			ArrayList<FeatureWeightPair> child2 = deepCopyIndividual(tournamentSelection(tournamentSize).individual);

			if (isLogging)
				log("mating " + getIndividualAsStr(child1) + " with " + getIndividualAsStr(child2));

			uniformCrossover(child1, child2);

			if (isLogging)
				log("After crossover: child1 = " + getIndividualAsStr(child1) + ", child2 = "
						+ getIndividualAsStr(child2));

			mutate(child1);

			mutate(child2);
			if (isLogging)
				log("After mutation: child1 = " + getIndividualAsStr(child1) + ", child2 = "
						+ getIndividualAsStr(child2));

			// normalize(child1);
			// normalize(child2);
			// if (isLogging)
			// log("After normalization: child1 = " + getIndividualAsStr(child1)
			// + ", child2 = " + getIndividualAsStr(child2));

			children.add(child1);
			children.add(child2);

		}

		return children;
	}

	// does a deep copy of an individual's weights
	ArrayList<FeatureWeightPair> deepCopyIndividual(ArrayList<FeatureWeightPair> individual) {
		ArrayList<FeatureWeightPair> copy = new ArrayList<FeatureWeightPair>(individual.size());

		for (int i = 0; i < individual.size(); i++)
			copy.add(new FeatureWeightPair(individual.get(i).feature, individual.get(i).weight,
					individual.get(i).increasesHappiness));

		return copy;
	}

	// this method is unusued
	// returns information about the lowest, average and highest score on an
	// individual after playing numGames games, each game with random piece
	// sequences
	private FitnessAssessment assessFitness(ArrayList<FeatureWeightPair> individual) {
		ArrayList<Integer> scores = new ArrayList<>(numGames);
		ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(numGames);

		for (int game = 0; game < numGames; game++) {
			PlayerSkeleton player = new PlayerSkeleton();
			player.setFeatureWeightPairs(individual);
			player.setSeed(rng.nextLong()); // this is the actual seed used in
											// the sequence
			tasks.add(service.submit(player));
		}

		for (int game = 0; game < numGames; game++)
			try {
				scores.add(tasks.get(game).get());
			} catch (Exception e) {
				log("the future computing a game's score was interupted");
				System.out.println("the future computing a game's score was interupted");
			}

		int lowest = scores.get(0);
		int highest = lowest;
		int sum = 0;

		for (Integer score : scores) {
			if (score > highest) {
				highest = score;
			}
			if (score < lowest) {
				lowest = score;
			}
			sum += score;
		}

		return new FitnessAssessment(individual, lowest, sum / numGames, highest);
	}

	// returns information about the lowest, average and highest score on the
	// entire population after playing numGames games, each game with random
	// piece sequences
	private void assessFitnessOfPopulation() {
		fitnessResults.clear();

		ArrayList<Integer> scores = new ArrayList<>(numGames * population.size());
		ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(numGames * population.size());
		long[] seeds = new long[numGames];

		for (int game = 0; game < numGames; game++)
			seeds[game] = rng.nextLong();

		for (ArrayList<FeatureWeightPair> individual : population) {
			for (int game = 0; game < numGames; game++) {
				PlayerSkeleton player = new PlayerSkeleton();
				player.setFeatureWeightPairs(individual);
				player.setSeed(seeds[game]);
				tasks.add(service.submit(player));
			}
		}

		for (int i = 0; i < population.size(); i++) {
			for (int game = 0; game < numGames; game++)
				try {
					int index = (i * numGames) + game;
					scores.add(tasks.get(index).get());
				} catch (Exception e) {
					log("the future computing a game's score was interupted");
					System.out.println("the future computing a game's score was interupted");
				}
		}

		for (int individual = 0; individual < population.size(); individual++) {
			int startIndex = individual * numGames; // the individual's scores
													// are in indices
													// [startIndex,
													// startIndex+numGames)
			int lowest = scores.get(startIndex); // the first game played by
													// this individual
			int highest = lowest;
			int sum = lowest;
			// check second game and up
			for (int game = 1; game < numGames; game++) {
				int score = scores.get(startIndex + game);
				
				if (score > highest) {
					highest = score;
				}
				if (score < lowest) {
					lowest = score;
				}
				sum += score;
			}

			fitnessResults.add(new FitnessAssessment(population.get(individual), lowest, sum / numGames, highest));
		}

	}

	// returns the fitness assessment of the individual being selected through
	// tournament selection
	private FitnessAssessment tournamentSelection(int tournamentSize) {
		int best = randomInt(0, population.size());
		FitnessAssessment bestFitness = fitnessResults.get(best);

		for (int i = 2; i <= tournamentSize; i++) {
			int next = randomInt(0, population.size());
			// System.out.println("i = " + i + " best = " + best + " next = " +
			// next);

			// ensure that the next individual selected is different from the
			// existing best
			while (best == next)
				next = randomInt(0, fitnessResults.size());

			FitnessAssessment fitness = fitnessResults.get(next);

			if (fitness.compareTo(bestFitness) > 0) {
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
			n = nextGaussian(0.0f, mutationSigma);
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
			// error opening the log file - just get rid of logging so it won't
			// print to the console while the user is running the program
			LogManager.getLogManager().reset();
		}
	}

	private static void log(String msg) {
		if (isLogging)
			logger.log(Level.INFO, msg);
	}

	public static void writeToCSV(FitnessAssessment fitness) {
		try {
			for (FeatureWeightPair feature : fitness.individual) {
				writer.append(String.valueOf(feature.weight));
				writer.append(',');
			}
			writer.append(String.valueOf(fitness.lowest));
			writer.append(',');
			writer.append(String.valueOf(fitness.average));
			writer.append(',');
			writer.append(String.valueOf(fitness.highest));
			writer.append('\n');
			writer.flush();
		} catch (Exception e) {
			// TODO: handle exception
		}
	}

	public static void main(String[] args) {
		int cores = Runtime.getRuntime().availableProcessors();
		System.out.println("Running genetic algorithm on a machine with " + cores + " logical cores:");
		Scanner sc = new Scanner(System.in);
		String filename;
		System.out.print("\nEnter log filename: ");
		filename = sc.nextLine();
		loggerInit(filename);
		int elites = 0, games = 0, populationSize = 0, tournamentSize = 0, generations = 0, convergenceThreshhold = 0,
				threads = 0, repeats = 0;
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
		System.out.print(
				"\nNumber of threads to use (aid for collecting report data. should not be > number of cores): ");
		threads = sc.nextInt();
		System.out.print(
				"\nEnter number of times to repeat training with above settings (for obtaining training times for report mostly): ");
		repeats = sc.nextInt();
		System.out.println("running " + repeats + " times:");

		long totalTime = 0;

		for (int i = 0; i < repeats; i++) {
			try {
				writer = new FileWriter("weightsData.csv");
				System.out.println("csv file created");
				writer.append("MaxHeight");
				writer.append(',');
				writer.append("MeanHeight");
				writer.append(',');
				writer.append("NumHoles");
				writer.append(',');
				writer.append("RowTransitions");
				writer.append(',');
				writer.append("StdDevHeight");
				writer.append(',');
				writer.append("SumOfPitDepths");
				writer.append(',');
				writer.append("WellSums");
				writer.append(',');
				writer.append("LowestScore");
				writer.append(',');
				writer.append("AverageScore");
				writer.append(',');
				writer.append("HighestScore");
				writer.append('\n');
				GeneticAlgorithm ga = new GeneticAlgorithm(crossoverRate, elites, games, mutationSigma, populationSize,
						tournamentSize, threads);
				long startTime = System.currentTimeMillis();
				FitnessAssessment result = ga.trainFor(generations, convergenceThreshhold);
				long endTime = System.currentTimeMillis();
				long timeElapsed = endTime - startTime;
				totalTime += timeElapsed;
				System.out.println("Training complete. The best individual is ");
				System.out.println(result);
				System.out.println("GA took: " + timeElapsed + " ms (" + timeElapsed / 1000 + " seconds)");

				writer.flush();
				writer.close();
			} catch (IOException e) {

			}
		}

		System.out.println("On average, training with the above settings take " + totalTime / repeats + " ms");
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
*/

/* particle.java

//each particle represents a candidate solution

public class Particle {

	public float[] position; // the current particle's position
	public float[] velocity;
	public float[] personalBestVector; // the position/weight vector where it
										// achieved the best score
	public ParticleSwarmOptimization.FitnessAssessment personalBestFitness;

	// the position vector supplied is used to determine the number of
	// features/coordinates
	// initializes a particle with an initial position and the 0 velocity vector
	public Particle(float[] position) {
		this.position = position;
		velocity = new float[position.length];
		personalBestVector = new float[position.length];

		for (int i = 0; i < position.length; i++) {
			velocity[i] = 0.0f;
			personalBestVector[i] = position[i];
		}

		personalBestFitness = new ParticleSwarmOptimization.FitnessAssessment(position, 0, 0, 0);
	}

	// initializes particle with position and velocity
	public Particle(float[] position, float[] velocity) {
		this.position = position;
		this.velocity = velocity;
		personalBestVector = new float[position.length];

		for (int i = 0; i < position.length; i++)
			personalBestVector[i] = position[i];

		personalBestFitness = new ParticleSwarmOptimization.FitnessAssessment(position, 0, 0, 0);
	}

	public void setBestFitnessVector(float[] pos) {
		for (int i = 0; i < pos.length; i++)
			personalBestVector[i] = pos[i];
	}

	// returns a deep copy of position
	public float[] getCopyOfPosition() {
		float[] posCopy = new float[position.length];

		for (int i = 0; i < position.length; i++)
			posCopy[i] = position[i];

		return posCopy;
	}
}
*/

/* ParticleSwarmOptimization.java

//main reference: essentials of metaheuristics, page 59
//this implementation is essentially straight from the book minus informants and particle jump size
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
	// since particles are just vectors in r^n, the correspondence between the
	// i^{th} coordinate
	// and the corresponding features and its min and max allowed weights is
	// maintained here
	private ArrayList<FeatureBoundaryPair> features;
	private ArrayList<Particle> population;
	private int populationSize; // size of the swarm
	// inertia is the proportion of velocity to be retained; or its "momentum"
	// so that the particle's movement isn't too erratic
	// Low settings favour exploitation (local search) while a high value
	// favours global search.
	// A decaying inertia weight can be used, with the aim of favoring global
	// search at the start of the algorithm and local search later.
	// If inertia is not reduced with time, choose a value in [0.8, 1.2]
	// see http://tracer.uc3m.es/tws/pso/parameters.html
	private float startInertia, endInertia, inertia;
	// the cognitive component is the proportion of the particle's personal best
	// to be retained
	// allows a particle's memory of where it has previously found good
	// solutions to influence its movement direction
	// If this is large, particles tend to move more towards their own personal
	// bests rather than towards global bests.
	// This breaks the swarm into a lot of separate hill-climbers rather than a
	// joint searcher.
	private float cognitive;
	// the social component is the proportion of global best to be retained.
	// allows the knowledge of other swarm members (i.e. where other members of
	// the swarm have found good solutions) to influence a particle's movement
	// If this is large, particles tend to move more towards the best known
	// region. This converts the algorithm into one large hill-climber
	// the book suggests setting it to 0, other sources disagree on this; we'll
	// have to experiment
	private float social;
	// other sources on parameter selection:
	// http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=4019136&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D4019136
	// http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=5957817&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D5957817
	// http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=870279&tag=1&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D870279%26tag%3D1

	public int numGames = 0; // number of games to run to assess fitness of an
								// individual
	private Random rng;
	private ArrayList<FitnessAssessment> fitnessResults;
	private ExecutorService service = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

	ParticleSwarmOptimization(int populationSize, int games, float cognitive, float startInertia, float endInertia,
			float social) {
		rng = new Random();
		this.populationSize = populationSize;
		this.numGames = games;
		this.cognitive = cognitive;
		this.startInertia = startInertia;
		this.endInertia = endInertia;
		inertia = startInertia;
		this.social = social;
		population = new ArrayList<Particle>(populationSize);
		fitnessResults = new ArrayList<FitnessAssessment>(populationSize);
		features = new ArrayList<FeatureBoundaryPair>();
		// add features here
		// features.add(new FeatureBoundaryPair(new PlayerSkeleton.Bumpiness(),
		// -10.0f, 0.0f));
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
		for (int i = 0; i < population.size(); i++) {
			str += i + ": " + getIndividualAsStr(i);
			str += "\r\n";
		}
		log(str);
	}

	private String getIndividualAsStr(Particle individual) {
		String str = "";

		for (int i = 0; i < individual.position.length; i++)
			str += individual.position[i] + ", ";

		str += " velocity: ";

		for (int i = 0; i < individual.position.length; i++)
			str += individual.velocity[i] + ", ";

		str += " personal best: ";

		for (int i = 0; i < individual.position.length; i++)
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
		for (int i = 0; i < populationSize; i++)
			population.add(generateRandomStationaryParticle());
	}

	// generates a random particle with a 0 velocity vector
	private Particle generateRandomStationaryParticle() {
		float[] pos = new float[features.size()];

		for (int i = 0; i < features.size(); i++)
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

		for (int i = 0; i < vec.length; i++)
			vec[i] /= length;
	}

	// run the algorithm for a specified number of generations
	public FitnessAssessment trainFor(int generations, int convergenceThreshhold) {
		System.out.println("Training for " + generations + " generations with convergence threshhold of "
				+ convergenceThreshhold + " generations:");
		log("training for " + generations + " generations with convergence threshhold of " + convergenceThreshhold
				+ " generations:");

		float[] bestPosition = population.get(0).getCopyOfPosition();
		// initialize to the worst possible score
		FitnessAssessment bestFitness = new FitnessAssessment(bestPosition, 0, 0, 0);
		int generationBestWasFound = 0;
		boolean convergence = false;
		float changeInInertia = (endInertia - startInertia) / generations;

		for (int generation = 0; generation < generations && convergence == false; generation++) {
			System.out.println("currently on generation " + generation);
			long[] seeds = new long[numGames];

			for (int game = 0; game < numGames; game++)
				seeds[game] = rng.nextLong();

			bestFitness = assessFitness(seeds, bestPosition);
			assessFitnessOfPopulation(seeds);

			for (int i = 0; i < population.size(); i++) {
				Particle particle = population.get(i);
				FitnessAssessment curFitness = fitnessResults.get(i);
				if (curFitness.compareTo(bestFitness) > 0) {
					bestFitness = curFitness;
					bestPosition = particle.getCopyOfPosition();
					generationBestWasFound = generation;
				}

				if (curFitness.compareTo(particle.personalBestFitness) > 0) {
					particle.personalBestFitness = curFitness;
					particle.setBestFitnessVector(particle.position); // we need
																		// a
																		// deep
																		// copy
																		// so
																		// position
																		// changes
																		// don't
																		// affect
																		// this
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
			System.out.println("Training stopped because convergence was detected with no improvement after "
					+ convergenceThreshhold + " generations");
			log("Training stopped because convergence was detected with no improvement after " + convergenceThreshhold
					+ " generations");
		}
		service.shutdown();

		return bestFitness;
	}

	// assess fitness of an individual
	private FitnessAssessment assessFitness(long[] seeds, float[] pos) {
		ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(numGames);
		ArrayList<Integer> scores = new ArrayList<Integer>(numGames);

		for (int game = 0; game < numGames; game++) {
			PlayerSkeleton player = getNewPlayerUsingWeights(pos);
			player.setSeed(seeds[game]);
			tasks.add(service.submit(player));
		}

		// retrieve scores for best vector
		for (int game = 0; game < numGames; game++)
			try {
				scores.add(tasks.get(game).get());
			} catch (Exception e) {
				log("the future computing a game's score was interupted");
				System.out.println("the future computing a game's score was interupted");
			}

		int lowest = scores.get(0); // the first game played by best vector
		int highest = lowest;
		int sum = lowest;
		// check second game and up
		for (int game = 1; game < numGames; game++) {
			int score = scores.get(game);
			if (score > highest) {
				highest = score;
			}
			if (score < lowest) {
				lowest = score;
			}
			sum += score;
		}

		return new FitnessAssessment(pos, lowest, sum / numGames, highest);
	}

	// returns a PlayerSkeleton initialized with a weight vector
	private PlayerSkeleton getNewPlayerUsingWeights(float[] weights) {
		PlayerSkeleton player = new PlayerSkeleton();
		int numFeatures = features.size();
		ArrayList<FeatureWeightPair> fwPairs = new ArrayList<FeatureWeightPair>(numFeatures);

		for (int i = 0; i < numFeatures; i++)
			fwPairs.add(new FeatureWeightPair(features.get(i).feature, weights[i], false)); // the
																							// third
																							// argument
																							// doesn't
																							// really
																							// matter
																							// here;
																							// only
																							// used
																							// by
																							// GA

		player.setFeatureWeightPairs(fwPairs);
		return player;
	}

	// computes the fitness of all particles
	private void assessFitnessOfPopulation(long[] seeds) {
		fitnessResults.clear();

		ArrayList<Integer> scores = new ArrayList<>(numGames * population.size());
		ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(numGames * population.size());

		for (Particle individual : population) {
			for (int game = 0; game < numGames; game++) {
				PlayerSkeleton player = getNewPlayerUsingWeights(individual.position);
				player.setSeed(seeds[game]);
				tasks.add(service.submit(player));
			}
		}

		for (int i = 0; i < population.size(); i++) {
			for (int game = 0; game < numGames; game++)
				try {
					int index = (i * numGames) + game;
					scores.add(tasks.get(index).get());
				} catch (Exception e) {
					log("the future computing a game's score was interupted");
					System.out.println("the future computing a game's score was interupted");
				}
		}

		for (int individual = 0; individual < population.size(); individual++) {
			int startIndex = individual * numGames; // the individual's scores
													// are in indices
													// [startIndex,
													// startIndex+numGames)
			int lowest = scores.get(startIndex); // the first game played by
													// this individual
			int highest = lowest;
			int sum = lowest;
			// check second game and up
			for (int game = 1; game < numGames; game++) {
				int score = scores.get(startIndex + game);
				if (score > highest) {
					highest = score;
				}
				if (score < lowest) {
					lowest = score;
				}
				sum += score;
			}

			fitnessResults
					.add(new FitnessAssessment(population.get(individual).position, lowest, sum / numGames, highest));
		}

		tasks.clear();
		scores.clear();
		// update the score of the personal best positions of the particle
		// we can't reuse the old value as the personal best would perform
		// differently across generations
		for (Particle individual : population) {
			for (int game = 0; game < numGames; game++) {
				PlayerSkeleton player = getNewPlayerUsingWeights(individual.personalBestVector);
				player.setSeed(seeds[game]);
				tasks.add(service.submit(player));
			}
		}

		for (int i = 0; i < population.size(); i++) {
			for (int game = 0; game < numGames; game++)
				try {
					int index = (i * numGames) + game;
					scores.add(tasks.get(index).get());
				} catch (Exception e) {
					log("the future computing a game's score was interupted");
					System.out.println("the future computing a game's score was interupted");
				}
		}

		for (int individual = 0; individual < population.size(); individual++) {
			int startIndex = individual * numGames; // the individual's scores
													// are in indices
													// [startIndex,
													// startIndex+numGames)
			int lowest = scores.get(startIndex); // the first game played by
													// this individual
			int highest = lowest;
			int sum = lowest;
			// check second game and up
			for (int game = 1; game < numGames; game++) {
				int score = scores.get(startIndex + game);
				if (score > highest) {
					highest = score;
				}
				if (score < lowest) {
					lowest = score;
				}

				sum += score;
			}

			FitnessAssessment f = new FitnessAssessment(population.get(individual).personalBestVector, lowest,
					sum / numGames, highest);
			population.get(individual).personalBestFitness = f;
		}

	}

	// determines the velocity of particles in this generation
	// bestSoFar is the coordinates of best particle found so far
	private void determineVelocity(float[] bestSoFar) {
		for (Particle p : population)
			for (int i = 0; i < features.size(); i++) {
				float c1 = randomFloat(0.0f, cognitive);
				float c2 = randomFloat(0.0f, social);
				float diffFromPersonalBest = p.personalBestVector[i] - p.position[i];
				float diffFromBestFoundSoFar = bestSoFar[i] - p.position[i];
				p.velocity[i] = (inertia * p.velocity[i]) + (c1 * diffFromPersonalBest) + (c2 * diffFromBestFoundSoFar);

			}
	}

	// mutates the particles by applying their velocities to their current
	// position
	private void mutate() {
		for (Particle p : population)
			for (int i = 0; i < features.size(); i++) {
				float newCoord = p.position[i] + p.velocity[i];
				newCoord = clamp(newCoord, features.get(i).min, features.get(i).max); // prevent
																						// particles
																						// from
																						// escaping
																						// boundaries
				p.position[i] = newCoord;
			}
	}

	// clamps a value if it is outside [min, max] so that it falls inside this
	// range
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
			// error opening the log file - just get rid of logging so it won't
			// print to the console while the user is running the program
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
		ParticleSwarmOptimization pso = new ParticleSwarmOptimization(populationSize, games, cognitive, startInertia,
				endInertia, social);
		long startTime = System.currentTimeMillis();
		FitnessAssessment result = pso.trainFor(generations, convergenceThreshhold);
		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("Training complete. The best individual is ");
		System.out.println(result);
		System.out.println("PSO took: " + totalTime + " ms (" + totalTime / 1000 + " seconds)");
		sc.close();
	}

	// a utility class that associates a feature with the minimum and maximum
	// weight
	// this constrains the search space
	private class FeatureBoundaryPair {
		float min = 0.0f, max = 0.0f;
		FeatureFunction feature;

		FeatureBoundaryPair(FeatureFunction f, float min, float max) {
			assert (min <= max);
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
*/

