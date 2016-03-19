import java.util.*;

public class PlayerSkeleton {

	public static final int COLS = 10;
	public static final int ROWS = 21;
	public static final int N_PIECES = 7;

	//implement this function to have a working system
	public int pickMove(State s, int[][] legalMoves) {
		// state is current state
		// legalMoves is legalMoves for next piece, 2D array [numLegalMoves][0 = orient/ 1 = slot]
		// for each legalMove, apply it to state and get the next state (see what each new state looks like)
		// for each potential new state, use the GA to get the ratings to determine best move
		// choose the best move and return it

		
		System.out.println(s.getNextPiece());

		Random rand = new Random(); // just for fun, better than the original given by prof. Comment away to see the original simulation given by prof
		int randMove = rand.nextInt(legalMoves.length); // just for fun, better than the original given by prof. Comment away to see the original simulation given by prof

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
		
		return randMove; // return moveNumber in legalMoves[moveNum][orient/slot]
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

	// function to return agg height: sum of height of each column
	private static int aggHeight(State s) {
		int aggHeight = 0;

		int[] top = s.getTop();

		for (int i = 0; i < COLS; i++) {
			aggHeight += top[i];
		}

		return aggHeight;
	}

	// function to return numCompleteLines: unsure
	private static int numCompleteLines(State s) {
		return s.getRowsCleared(); // this gets the cumulative total number of rows cleared so far. Do we want this? 
		// Or do y'all want to calculate manually for each new state (not cumulative, only count rows cleared by current action)
		// i'm not sure if the State.java allow us to maintain a state with complete lines (for us to count) without actually executing it on the animation :/
	}

	// function to return numHoles: A hole is an empty space such that there is at least one tile in the same column above it.
	private static int numHoles(State s) {
		int numHoles = 0;

		int[][] field = s.getField();
		int[] top = s.getTop();

		for (int i = 0; i < COLS; i++) {
			for (int j = 0; j < top[i]; j++) {
				if (field[j][i] == 0) {
					numHoles++;
				}
			}
			// System.out.println("numHoles: " + numHoles);
		}

		return numHoles;
	}

	// function to return bumpiness: summing up the absolute differences between all two adjacent columns
	private static int bumpiness(State s) {
		int bumpiness = 0;

		int[] top = s.getTop();

		for (int i = 0; i < COLS - 1; i++) {
			bumpiness += Math.abs(top[i] - top[i+1]);
		}
		// System.out.println("bumpiness: " + bumpiness);

		return bumpiness;
	}

	// function to return max column height: DONE
	private static int maxColumnHeight(State s) {
		int maxColumnHeight = 0;

		int[] top = s.getTop();

		for (int i = 0; i < COLS; i++) {
			// System.out.println("top: " + top[i]);
			if (top[i] > maxColumnHeight) {
				maxColumnHeight = top[i];
			}
			// System.out.println("maxColumnHeight: " + maxColumnHeight);
		}
		return maxColumnHeight;
	}
	
}
