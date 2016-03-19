import java.util.*;
import java.util.*;

public class PlayerSkeleton {

	public static final int COLS = 10;
	public static final int ROWS = 21;
	public static final int N_PIECES = 7;
private ArrayList<FeatureFunction> features;
private ArrayList<Float> weights;

	public PlayerSkeleton() {
		features = new ArrayList<FeatureFunction>();
		weights = new ArrayList<Float>();
		// to standardize, lets add the features and their corresponding weights in alphabetical order?
		// these weights are just arbitrarily chosen starting points, we can tweak them later with the GA 
		features.add(new AggHeight());
		weights.add(new Float(-0.5));
		features.add(new Bumpiness());
		weights.add(new Float(-0.25));
		features.add(new MaxHeight());
		weights.add(new Float(-0.5));
features.add(new NumHoles());
weights.add(new Float(-0.5));
	}
	
	// returns  h(n) for the given state
	private float evaluate(ImprovedState s) {
		float sum = 0.0f;
		
		for (int i=0; i<features.size(); i++)
			sum += features.get(i).evaluate(s) * weights.get(i);
		
			return sum;
	}
	
	//implement this function to have a working system
	public int pickMove(State s, int[][] legalMoves) {
		// s is the current state
		// legalMoves is legalMoves for next piece, 2D array [numLegalMoves][0 = orient/ 1 = slot]
		
		ImprovedState currentState = new ImprovedState(s);

		int bestMove = 0;
		boolean isNonLosingMoveFound = false;
		float bestValue = 0.0f;
		
		// find the first legal move corresponding to a non-losing situation, and assume that to be the best move
		for (int move = 0; move < legalMoves.length && isNonLosingMoveFound == false; move++) {
			ImprovedState resultingState = currentState.tryMove(move);
			isNonLosingMoveFound = !resultingState.hasLost();

			if (isNonLosingMoveFound) {
				bestMove = move;
				bestValue = evaluate(resultingState);
			}
		}
			
		if (isNonLosingMoveFound == false)
			return 0; // if we'll die anyway, it doesn't matter which move we do
		
		// now see if we can find better moves
		for (int move = bestMove; move < legalMoves.length; move++) {
			ImprovedState resultingState = currentState.tryMove(move);
		float utility = evaluate(resultingState);	
		
		if (utility > bestValue) {
			bestValue = utility;
			bestMove = move;
		}
		}
		
		return bestMove;

		/* 
		System.out.println(s.getNextPiece());

		Random rand = new Random(); // just for fun, better than the original given by prof. Comment away to see the original simulation given by prof
		int randMove = rand.nextInt(legalMoves.length); // just for fun, better than the original given by prof. Comment away to see the original simulation given by prof
*/
		
		/* For debugging: to see what is stored in legalMoves

		for (int i = 0; i < legalMoves.length; i++) {
			for (int j = 0; j < legalMoves[i].length; j++) {
				System.out.println("legalMoves[" + i + "][" + j + "]: " + legalMoves[i][j]);
			}
		}

		System.out.println("legalMoves[] length: " + legalMoves.length);
		for (int i = 0; i < legalMoves.length; i++) {
			System.out.println("legalMoves[][] length: " + legalMoves[i].length);
		}
		*/
		
// 		return randMove; // return moveNumber in legalMoves[moveNum][orient/slot]
	}
	
	public static void main(String[] args) {
		State s = new State();
		new TFrame(s);
		PlayerSkeleton p = new PlayerSkeleton();
		while(!s.hasLost()) {
			s.makeMove(p.pickMove(s,s.legalMoves()));
			s.draw();
			s.drawNext(0,0);
			try {
				Thread.sleep(300);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		System.out.println("You have completed "+s.getRowsCleared()+" rows.");
	}

	// GA: to be implemented
	// Population: all the states generated from current state, given the legal moves (population should be passed in as an argument when GA is called from pickMove())
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

	// 	// function to return numCompleteLines: unsure
	private static int numCompleteLines(State s) {
		return s.getRowsCleared(); // this gets the cumulative total number of rows cleared so far. Do we want this? 
		// Or do y'all want to calculate manually for each new state (not cumulative, only count rows cleared by current action)
		// i'm not sure if the State.java allow us to maintain a state with complete lines (for us to count) without actually executing it on the animation :/
	}


	// returns aggregate height for all columns
			private class AggHeight implements FeatureFunction {
				@Override
		public float evaluate(ImprovedState s)  {
				int aggHeight = 0;

				int[] top = s.getTop();

				for (int i = 0; i < State.COLS; i++) {
					aggHeight += top[i];
				}

				return (float) aggHeight;
			}
			}
			
			// returns number of holes. A hole is an empty space such that there is at least one tile in the same column above it
	private class NumHoles implements FeatureFunction {
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
			// System.out.println("numHoles: " + numHoles);
			return numHoles;
		}
	}

	// calculates bumpiness, the sum of the absolute differences between heights of consecutive adjacent columns
	private class Bumpiness implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int bumpiness = 0;

			int[] top = s.getTop();

			for (int i = 0; i < State.COLS - 1; i++) {
				bumpiness += Math.abs(top[i] - top[i+1]);
			}
			// System.out.println("bumpiness: " + bumpiness);

			return bumpiness;
		}
	}

		// maximum column height heuristic: DONE
	private class MaxHeight implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int maxColumnHeight = 0;

			int[] top = s.getTop();

			for (int i = 0; i < State.COLS; i++) {
				// System.out.println("top: " + top[i]);
				if (top[i] > maxColumnHeight) {
					maxColumnHeight = top[i];
				}
				// System.out.println("maxColumnHeight: " + maxColumnHeight);
			}
			return maxColumnHeight;
		}
	}
		
}

