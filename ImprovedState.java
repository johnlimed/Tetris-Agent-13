import java.util.Random;

// this is a copy of the state file in the skeleton 
// with the ability to simulate the field when trying out different moves, and to set the seed used for producing sequences    
// usage: 
// construct an ImprovedState object with a State object so it has a copy of the current state
// then, use tryMove to get all the states resulting from trying out all possible moves  


public class ImprovedState {
	public static final int COLS = 10;
	public static final int ROWS = 21;
	public static final int N_PIECES = 7;
	// this rng is used to ensure identical piece sequences
private Random rng = new java.util.Random();

	public boolean lost = false, prevLost;
	
	private int turn = 0, prevTurn;
	private int cleared = 0, prevCleared; 
	
	//each square in the grid - int means empty - other values mean the turn it was placed
	private int[][] curField = new int[ROWS][COLS]; // the current state of the grid
	private int[][] prevField; // the previous field, for undo functionality
	
	//top row+1 of each column
	//0 means empty
	private int[] top = new int[COLS];
	private int[] prevTop;
	
	
	//number of next piece
	protected int nextPiece, prevPiece;
	
	
	//all legal moves - first index is piece type - then a list of 2-length arrays
	protected static int[][][] legalMoves = new int[N_PIECES][][];
	
	//indices for legalMoves
	public static final int ORIENT = 0;
	public static final int SLOT = 1;
	
	//possible orientations for a given piece type
	protected static int[] pOrients = {1,2,4,4,4,2,2};
	
	//the next several arrays define the piece vocabulary in detail
	//width of the pieces [piece ID][orientation]
	protected static int[][] pWidth = {
			{2},
			{1,4},
			{2,3,2,3},
			{2,3,2,3},
			{2,3,2,3},
			{3,2},
			{3,2}
	};
	//height of the pieces [piece ID][orientation]
	private static int[][] pHeight = {
			{2},
			{4,1},
			{3,2,3,2},
			{3,2,3,2},
			{3,2,3,2},
			{2,3},
			{2,3}
	};
	private static int[][][] pBottom = {
		{{0,0}},
		{{0},{0,0,0,0}},
		{{0,0},{0,1,1},{2,0},{0,0,0}},
		{{0,0},{0,0,0},{0,2},{1,1,0}},
		{{0,1},{1,0,1},{1,0},{0,0,0}},
		{{0,0,1},{1,0}},
		{{1,0,0},{0,1}}
	};
	private static int[][][] pTop = {
		{{2,2}},
		{{4},{1,1,1,1}},
		{{3,1},{2,2,2},{3,3},{1,1,2}},
		{{1,3},{2,1,1},{3,3},{2,2,2}},
		{{3,2},{2,2,2},{2,3},{1,2,1}},
		{{1,2,2},{3,2}},
		{{2,2,1},{2,3}}
	};
	
	//initialize legalMoves array for each piece
	{
		//for each piece type
		for(int i = 0; i < N_PIECES; i++) {
			//figure number of legal moves
			int n = 0;
			for(int j = 0; j < pOrients[i]; j++) {
				//number of locations in this orientation
				n += COLS+1-pWidth[i][j];
			}
			//allocate space
			legalMoves[i] = new int[n][2]; // n is number of legalMoves for this piece, second[] has 2 values containing: ORIENT/SLOT
			//for each orientation
			n = 0;
			for(int j = 0; j < pOrients[i]; j++) {
				//for each slot
				for(int k = 0; k < COLS+1-pWidth[i][j];k++) {
					legalMoves[i][n][ORIENT] = j; // i is the piece type, j = number of possible orients
					legalMoves[i][n][SLOT] = k;
					n++;
				}
			}
		}
	
	}
	
	public ImprovedState() {
		nextPiece = randomPiece();
		initPrevState();
	}
	
	// initialize all the variables that comprise the previous state
	private void initPrevState() {
		prevLost = lost;
		prevTurn = turn;
		prevField= new int[ROWS][COLS];
		
		for (int row = 0; row < ROWS; row++)
			System.arraycopy(curField[row], 0, prevField[row], 0, COLS);
		
		prevTop = new int[COLS];
		System.arraycopy(top, 0, prevTop, 0, COLS);
		
		prevPiece = nextPiece;
		prevCleared = cleared;
	}
	
	// does deep copy of the relevant state information
		public ImprovedState(State s) {
			this(s.getField(), s.getTop(), s.lost, s.getNextPiece(), s.getRowsCleared(), s.getTurnNumber(), true);
		}
		
				public ImprovedState(ImprovedState s) {
					this(s.getField(), s.getTop(), s.lost, s.getNextPiece(), s.getRowsCleared(), s.getTurnNumber(), true);
				}
				
		public void setSeed(long seed) {
			rng.setSeed(seed);
		}
		
	//random integer, returns 0-6
		private int randomPiece() {
			/*
			int rand = rng.nextInt((8 - 0) + 1) + 0;
			
			if (rand == 7)
				return 5;
			else if (rand == 8)
				return 6;
			else
				return rand;
			*/
			
			return (int)(rng.nextDouble()*N_PIECES);
		}

public int pickNextPiece() {
	//pick a new piece
	prevPiece = nextPiece;
	nextPiece = randomPiece();
return nextPiece;
}

public void setNextPiece(int piece) {
	prevPiece = nextPiece;
	nextPiece = piece;
}
	public int[][] getField() {
		return curField;
	}

	public int[] getTop() {
		return top;
	}

    public static int[] getpOrients() {
        return pOrients;
    }
    
    public static int[][] getpWidth() {
        return pWidth;
    }

    public static int[][] getpHeight() {
        return pHeight;
    }

    public static int[][][] getpBottom() {
        return pBottom;
    }

    public static int[][][] getpTop() {
        return pTop;
    }

	public int getNextPiece() {
		return nextPiece;
	}

	public boolean hasLost() {
		return lost;
	}
	
	public int getRowsCleared() {
		return cleared;
	}
	
	public int getTurnNumber() {
		return turn;
	}
	
	// the deepCopy parameter controls whether a deep copy of the objects are done, which is required for modification not to affect other objects
	// Being able to control whether copies are deep or not is primarily a performance optimization
	private ImprovedState(int[][] field, int[] top, boolean lost, int nextPiece, int cleared, int turn, boolean deepCopy) {
		if (deepCopy) {
		for (int row = 0; row < ROWS; row++) {
		System.arraycopy(field[row], 0, this.curField[row], 0, COLS);
		}
		
		System.arraycopy(top, 0, this.top, 0, COLS);
		} else {
			this.curField = field;
			this.top = top;
		}
		
		this.nextPiece = nextPiece;
		this.lost = lost;
		this.cleared = cleared;
		this.turn = turn;
		initPrevState();
	}
	
	//gives legal moves for 
		public int[][] legalMoves() {
			return legalMoves[nextPiece];
		}
		
public int[][][] getLegalMovesForAllPieces() {
	return legalMoves;
}
		
		//returns a new object representing the state of the game if a particular move was made
		// move should be an index in the legal moves list
				// the returned ImprovedState object ccan then be fed into feature functions
		public ImprovedState tryMove(int move) {
			return tryMove(legalMoves[nextPiece][move]);
		}
		
		//make a move based on an array of orient and slot
		public ImprovedState tryMove(int[] move) {
			return tryMove(move[ORIENT],move[SLOT]);
		}
		
		private ImprovedState tryMove(int orient, int slot) {
			int[][] newField = new int[ROWS][COLS]; 
			
			for (int row = 0; row < ROWS; row++) {
				System.arraycopy(curField[row], 0, newField[row], 0, COLS);
				}
				
			int[] newTop = new int[COLS];
				System.arraycopy(top, 0, newTop, 0, COLS);
				int newCleared = cleared;
				int newTurn = turn + 1;
				
			//height if the first column makes contact
			int height = newTop[slot]-pBottom[nextPiece][orient][0];
			//for each column beyond the first in the piece
			for(int c = 1; c < pWidth[nextPiece][orient];c++) {
				height = Math.max(height,newTop[slot+c]-pBottom[nextPiece][orient][c]);
			}
			
			//check if game ended
			if(height+pHeight[nextPiece][orient] >= ROWS) {
				return new ImprovedState(newField, newTop, true, nextPiece, newCleared, newTurn, false);
			}

			//for each column in the piece - fill in the appropriate blocks
			for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
				
				//from bottom to top of brick
				for(int h = height+pBottom[nextPiece][orient][i]; h < height+pTop[nextPiece][orient][i]; h++) {
					newField[h][i+slot] = newTurn;
				}
			}
			
			//adjust top
			for(int c = 0; c < pWidth[nextPiece][orient]; c++) {
				newTop[slot+c]=height+pTop[nextPiece][orient][c];
			}
			
			int rowsCleared = 0;
			
			//check for full rows - starting at the top
			for(int r = height+pHeight[nextPiece][orient]-1; r >= height; r--) {
				//check all columns in the row
				boolean full = true;
				for(int c = 0; c < COLS; c++) {
					if(newField[r][c] == 0) {
						full = false;
						break;
					}
				}
				//if the row was full - remove it and slide above stuff down
				if(full) {
					rowsCleared++;
					newCleared++;
					//for each column
					for(int c = 0; c < COLS; c++) {

						//slide down all bricks
						for(int i = r; i < newTop[c]; i++) {
							newField[i][c] = newField[i+1][c];
						}
						//lower the top
						newTop[c]--;
						while(newTop[c]>=1 && newField[newTop[c]-1][c]==0)	newTop[c]--;
					}
				}
			}

			return new ImprovedState(newField, newTop, false, nextPiece, newCleared, newTurn, false);
		}

		// these functions execute the move and modifies the board; they're essentially copied directly from State
		//make a move based on the move index - its order in the legalMoves list
		public void makeMove(int move) {
			makeMove(legalMoves[nextPiece][move]);
		}
		
		//make a move based on an array of orient and slot
		public void makeMove(int[] move) {
			makeMove(move[ORIENT],move[SLOT]);
		}
		
		//returns false if you lose - true otherwise
		private boolean makeMove(int orient, int slot) {
			// make the current state the previous state
			for (int row = 0; row < ROWS; row++) {
				System.arraycopy(curField[row], 0, prevField[row], 0, COLS);
				}
			
			prevCleared = cleared;
			prevLost = lost;
			prevTurn = turn;
			prevPiece = nextPiece;
			turn++;
			System.arraycopy(top, 0, prevTop, 0, COLS);
			
			//height if the first column makes contact
			int height = top[slot]-pBottom[nextPiece][orient][0];
			//for each column beyond the first in the piece
			for(int c = 1; c < pWidth[nextPiece][orient];c++) {
				height = Math.max(height,top[slot+c]-pBottom[nextPiece][orient][c]);
			}
			
			//check if game ended
			if(height+pHeight[nextPiece][orient] >= ROWS) {
				lost = true;
				return false;
			}

			//for each column in the piece - fill in the appropriate blocks
			for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
				
				//from bottom to top of brick
				for(int h = height+pBottom[nextPiece][orient][i]; h < height+pTop[nextPiece][orient][i]; h++) {
					curField[h][i+slot] = turn;
				}
			}
			
			//adjust top
			for(int c = 0; c < pWidth[nextPiece][orient]; c++) {
				top[slot+c]=height+pTop[nextPiece][orient][c];
			}
			
			int rowsCleared = 0;
			
			//check for full rows - starting at the top
			for(int r = height+pHeight[nextPiece][orient]-1; r >= height; r--) {
				//check all columns in the row
				boolean full = true;
				for(int c = 0; c < COLS; c++) {
					if(curField[r][c] == 0) {
						full = false;
						break;
					}
				}
				//if the row was full - remove it and slide above stuff down
				if(full) {
					rowsCleared++;
					cleared++;
					//for each column
					for(int c = 0; c < COLS; c++) {

						//slide down all bricks
						for(int i = r; i < top[c]; i++) {
							curField[i][c] = curField[i+1][c];
						}
						//lower the top
						top[c]--;
						while(top[c]>=1 && curField[top[c]-1][c]==0)	top[c]--;
					}
				}
			}
	
			
			return true;
		}

		// undoes the effects of making the most recent move
		// only 1 level of undo is supported 
		public void undo() {
			// make the previous state the current
			int[][] tempField = curField;
			curField = prevField;
			prevField = tempField;

			int[] tempTop = top;
			top = prevTop;
			prevTop = tempTop;

			/*
			// swap the contents of prev and cur
			int[][] tempArr  = new int[ROWS][COLS];
					for (int row = 0; row < ROWS; row++) {
						System.arraycopy(curField[row], 0, tempArr[row], 0, COLS); // temp = cur
						System.arraycopy(prevField[row], 0, curField[row], 0, COLS); // cur = prev
						System.arraycopy(tempArr[row], 0, prevField[row], 0, COLS); // prev = temp
					}
					
						int[] tempTop = new int[COLS];
															System.arraycopy(top, 0, tempTop, 0, COLS); // temp = cur
															System.arraycopy(prevTop, 0, top, 0, COLS); // cur = prev
															System.arraycopy(tempTop, 0, prevTop, 0, COLS); // prev = temp
		*/
			
			int tempCleared = cleared;
			cleared = prevCleared;
						prevCleared = tempCleared;
			
						boolean tempLost = lost;
						lost = prevLost;
									prevLost = tempLost;
						
									int tempTurn = turn;
									turn = prevTurn;
												prevTurn = tempTurn;
												
												int tempPiece = nextPiece;
												nextPiece = prevPiece;
															prevPiece = tempPiece;
															
																												
		}

		// tests if the current state information in 2 ImprovedState objects are equal
		// This is mainly a debugging aid
		public boolean isCurStateEqual(ImprovedState other) {
boolean isCurStateEqual = lost == other.lost && turn == other.turn && cleared == other.cleared;
 
if (isCurStateEqual == false)
	return false;

for (int r=0; r<ROWS; r++)
	for (int c=0; c<COLS; c++)
		if (curField[r][c] != other.curField[r][c])
			return false;

for (int c=0; c<COLS; c++)
	if (top[c] != other.top[c])
		return false;

return true;
		}
		
				}