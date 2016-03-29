import java.util.*;


public class PlayerSkeleton {

	public static final int COLS = 10;
	public static final int ROWS = 21;
	public static final int N_PIECES = 7;

	private static long MAX_THREADS;
	private static long RUNNING_THREADS=0;
	private ArrayList<FeatureWeightPair> features;


	public PlayerSkeleton() {
		features = new ArrayList<FeatureWeightPair>();
	}


	public void setFeatureWeightPairs(ArrayList<FeatureWeightPair> features) {
		this.features = features;
	}

	// returns  h(n) for the given state
	private float evaluate(ImprovedState s) {
		float sum = 0.0f;

		for (int i=0; i<features.size(); i++) {
			FeatureWeightPair f = features.get(i);
			sum += f.feature.evaluate(s) * f.weight;
		}

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

		if (RUNNING_THREADS<MAX_THREADS) {
            long freeThreads = MAX_THREADS - RUNNING_THREADS;
            for(int i=0; i<freeThreads; i++) {
                Slave freeSlave = new Slave();
                freeSlave.start();
                RUNNING_THREADS++;
            }
        }


		// now see if we can find better moves
		for (int move = bestMove; move < legalMoves.length; move++) {
            ImprovedState resultingState = currentState.tryMove(move);
			if (resultingState.hasLost() == false) {
				float utility = evaluate(resultingState);

				if (utility > bestValue) {
					bestValue = utility;
					bestMove = move;
				}
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

	// plays a game , returning the number of rows completed
	// use the setFeatureWeightPairs function first
	public int playGame(boolean draw) {
		assert (features.isEmpty() == false); // must set some features to use first
		State s = new State();

		if (draw)
			new TFrame(s);

		while(!s.hasLost()) {
			s.makeMove(pickMove(s,s.legalMoves()));

			if (draw) {
				s.draw();
				s.drawNext(0,0);

				try {
					Thread.sleep(100);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}

		}

		return s.getRowsCleared();
	}

	public static void main(String[] args) {
        MAX_THREADS = Runtime.getRuntime().availableProcessors();
        System.out.println("number of processors: " + MAX_THREADS);
		PlayerSkeleton p = new PlayerSkeleton();
		// these weights are from a test run with 2 generations of the GA
		ArrayList<FeatureWeightPair> fwPairs = new ArrayList<FeatureWeightPair>();
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.AggHeight(), -1.5823991f));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.Bumpiness(), -0.33376628f));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.MaxHeight(), -0.19359511f));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.NumHoles(), -0.23551378f));
		p.setFeatureWeightPairs(fwPairs);
		System.out.println("You have completed "+p.playGame(true) +" rows.");
		/* 
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
		*/

	}


	// function to return numCompleteLines: unsure
	private static int numCompleteLines(State s) {
		return s.getRowsCleared(); // this gets the cumulative total number of rows cleared so far. Do we want this? 
		// Or do y'all want to calculate manually for each new state (not cumulative, only count rows cleared by current action)
		// i'm not sure if the State.java allow us to maintain a state with complete lines (for us to count) without actually executing it on the animation :/
	}


	// returns aggregate height for all columns
	public static class AggHeight implements FeatureFunction {
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
			// System.out.println("numHoles: " + numHoles);
			return numHoles;
		}
	}

	// calculates bumpiness, the sum of the absolute differences between heights of consecutive adjacent columns
	public static class Bumpiness implements FeatureFunction {
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
	public static class MaxHeight implements FeatureFunction {
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

