import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class PlayerSkeleton implements Callable<Integer> {

	public static final int COLS = 10;
	public static final int ROWS = 21;
	public static final int N_PIECES = 7;
	// multithreading of moves is currently not working; it returns different results from the sequential version, so uncommenting for now
    // private ExecutorService pickMoveThreadPool = Executors.newWorkStealingPool();
	private ArrayList<FeatureWeightPair> features;
private Long seed;

	public PlayerSkeleton() {
		features = new ArrayList<FeatureWeightPair>();
	}

	public void setFeatureWeightPairs(ArrayList<FeatureWeightPair> features) {
		this.features = features;
	}

	// returns  h(n) for the given state
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


	//implement this function to have a working system
	public int pickMove(State s, int[][] legalMoves) {
		int bestMove = 0;	
        float bestValue = Float.NEGATIVE_INFINITY;;
        /* 
if (isParallel) {
        List<Future<Float>> results = new ArrayList<Future<Float>>();
        Float[] value = new Float[legalMoves.length];
        ImprovedState currentState = new ImprovedState(s);
        for (int move=0; move<legalMoves.length; move++) {
            Slave slavePickMove = new Slave(currentState, move, features);
            results.add(pickMoveThreadPool.submit(slavePickMove));
        }
        
        int counter=0;
        for(Future<Float> result: results) {
            try {
                value[counter] = result.get();
            } catch (Exception e) {
                System.out.println("Error caught");
                pickMoveThreadPool.shutdown();
            }
            counter++;
        }
        bestValue = value[0];
        for (int i=0; i<value.length; i++) {
            // System.out.println("value: " + value[i] + " move: " + i);
            if (value[i] > bestValue) {
                bestValue = value[i];
                bestMove = i;
            }
        }
        
}
*/
        // non-parallel
ImprovedState currentState = new ImprovedState(s);
for (int move = 0; move < legalMoves.length; move++) {
	currentState.makeMove(move);
            float utility = evaluate(currentState);
currentState.undo();
            if (utility > bestValue) {
                bestValue = utility;
                bestMove = move;
            }
        }

return bestMove;
	}
 
	// picks the best move when running the game using ImprovedState
	public int pickMoveForImprovedState(ImprovedState s, int[][] legalMoves) throws InterruptedException {
		int bestMove = 0;	
        float bestValue = Float.NEGATIVE_INFINITY;;
        /*
int seq, par;
 
        {
        	ArrayList<Callable<Float>> tasks = new ArrayList<Callable<Float>>();
        List<Future<Float>> results;
        Float[] value = new Float[legalMoves.length];
        ImprovedState currentState = new ImprovedState(s);
        for (int move=0; move<legalMoves.length; move++) {
            Slave slavePickMove = new Slave(currentState, move, features);
            tasks.add(slavePickMove);
            // results.add(pickMoveThreadPool.submit(slavePickMove));
        }
        results = pickMoveThreadPool.invokeAll(tasks);
        int counter=0;
        for(Future<Float> result: results) {
            try {
                value[counter] = result.get();
            } catch (Exception e) {
                System.out.println("Error caught");
                pickMoveThreadPool.shutdown();
            }
            counter++;
        }
        
        bestValue = value[0];
        for (int i=0; i<value.length; i++) {
            // System.out.println("value: " + value[i] + " move: " + i);
            if (value[i] > bestValue) {
                bestValue = value[i];
                bestMove = i;
            }
        }
        par = bestMove;
}
		bestMove = 0;	
        bestValue = Float.NEGATIVE_INFINITY;;

        */
        
 // non-parallel implementation
// ImprovedState copy = new ImprovedState(s);
for (int move = 0; move < legalMoves.length; move++) {

	// ImprovedState resultingState = s.tryMove(move);
    // float utility = evaluate(resultingState);
// if (utility > bestValue) {
                // bestValue = utility;
                // bestMove = move;
            // }
           
    s.makeMove(move);
	float utility = evaluate(s);
	s.undo();
	
	// if (s.isCurStateEqual(copy) == false)
// 		System.out.println("current state not equals.");
	
	if (utility > bestValue) {
        bestValue = utility;
        bestMove = move;
    }
        }

// for checking if the results returned for both are different
// if (seq != par)
// System.out.println("Sequential implementation: " + seq + " parallel: " + par);

return bestMove;
	}

	// lookahead is the number of future pieces we will consider
	// for example, a value of 1 means that for the current state, execute each legal move m to get some new state. 
	// From this new state, calculate the highest utility each of the 7 pieces can have and average this out to get the utility for move m
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
		}
	 
	else {
		for (int move = 0; move < movesForAllPieces[nextPiece].length; move++) {
			ImprovedState future = s.tryMove(movesForAllPieces[nextPiece][move]);
			float sum = 0.0f; // sum of future utilities
			
			for (int piece=0; piece<ImprovedState.N_PIECES; piece++) {
				future.setNextPiece(piece);
				sum += pickMove(future, piece, lookahead-1).utility;
			}
			float avgFutureUtility = sum / ImprovedState.N_PIECES;
			if (avgFutureUtility > best) {
best = avgFutureUtility;
bestMove = move;
			}

		}
	}
		
	return new MoveUtilityPair(bestMove, best); 
	}
		
		// plays a game , returning the number of rows completed
	// use the setFeatureWeightPairs function first
	public int playGame(boolean alwaysDraw, boolean drawOnLoss) throws InterruptedException {
		assert (!features.isEmpty()); // must set some features to use first
		State s = new State();
		if (alwaysDraw)
			new TFrame(s);

		while(!s.hasLost()) {
			s.makeMove(pickMove(s,s.legalMoves()));
			if (alwaysDraw) {
				s.draw();
				s.drawNext(0,0);

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
			s.drawNext(0,0);
		}
		
		return s.getRowsCleared();
	}

	// sets the seed for playing a game with ImprovedState
	public void setSeed(long seed) {
		this.seed = seed;
	}
	
	// identical to the playGame function except that ImprovedState is used which is much more useful for training purposes
	// note that drawing isn't supported though with this one
		// use the setFeatureWeightPairs function first
		public int playGameWithImprovedState() throws InterruptedException {
			assert (!features.isEmpty()); // must set some features to use first
			ImprovedState s = new ImprovedState();
			if (seed != null) {
			s.setSeed(seed);
			}
			s.setSeed(1459523385737L);
			s.pickNextPiece();
						while(!s.hasLost()) {
			// for (int i=0; i<1; i++) {
				// s.makeMove(pickMoveForImprovedState(s,s.legalMoves()));
							s.makeMove(pickMove(s, s.getNextPiece(), 0).move);
				s.pickNextPiece();
				System.out.println(s.getRowsCleared());
								}
			
			return s.getRowsCleared();	
		}

	public static void main(String[] args) throws InterruptedException {
		PlayerSkeleton p = new PlayerSkeleton();
		// these weights are from a test run with 2 generations of the GA
		ArrayList<FeatureWeightPair> fwPairs = new ArrayList<FeatureWeightPair>();
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.AggHeight(), -0.4164334f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.Bumpiness(), -0.81970763f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.MaxHeight(), -0.93703204f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.NumHoles(), -5.9621563f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.MeanHeightDiff(), -4.838464f, false));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.NumRowsCleared(), 2.4920862f, true));
		fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.SumOfPitDepth(), -1.1674749f, false));
		p.setFeatureWeightPairs(fwPairs);
		
		long startTime = System.currentTimeMillis();
		System.out.println("You have completed "+p.playGameWithImprovedState() +" rows.");
		// System.out.println("You have completed "+p.playGame(false, true) +" rows.");
		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("PlayerSkeleton took: "+totalTime+"ms");
		
	}


		// Considered as a pit if the adjacent columns are >= 2. Depth = diff in height with the shortest adjacent col
	// should refactor this some more
	public static class SumOfPitDepth implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s)  {
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
				pitHeight = top[col+1];
				heightOfRightCol = top[col+2];

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
	public static class MeanHeightDiff implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s)  {
			int avgHeight = 0;
			int[] top = s.getTop();
			for (int i = 0; i < State.COLS; i++) {
				avgHeight += top[i];
			}
			return (float) avgHeight/State.COLS;
		}
    }
	
	public static class NumRowsCleared implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s)  {
			return (float) s.getRowsCleared();
		}
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

	// maximum column height heuristic
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
	// computes total row transitions. Row transitions happen when an empty cell is adjacent to a filled cell and vice versa
	public static class RowTransitions implements FeatureFunction {
		@Override
		public float evaluate(ImprovedState s) {
			int nRowTransitions = 0;
			int[][] field = s.getField();

			for (int r=0; r<State.ROWS; r++)
				for (int c=0; c<State.COLS-1; c++) {
					boolean isCurEmpty = field[r][c] == 0, isNextEmpty = field[r][c+1] == 0;
					
					if ((isCurEmpty && !isNextEmpty) || (!isCurEmpty && isNextEmpty))
						nRowTransitions++;
				}
			return nRowTransitions;
		}
	}
	
	// convenience class for associating a move index with its utility
private class MoveUtilityPair {
	public float utility = 0.0f;
	public int move = 0;
	
	MoveUtilityPair(int m, float u) {
		move = m;
		utility = u;
	}
}
}