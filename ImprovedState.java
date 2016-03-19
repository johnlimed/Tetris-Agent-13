// this is a copy of the state file in the skeleton 
// with the ability to simulate the field when trying out different moves    
// usage: 
// construct an ImprovedState object with a State object so it has a copy of the current state
// then, use tryMove to get all the states resulting from trying out all possible moves  
// finally, use the GA to pick the best one, then make the move on the actual State object

public class ImprovedState {
	public static final int COLS = 10;
	public static final int ROWS = 21;
	public static final int N_PIECES = 7;


	public boolean lost = false;
	
	private int turn = 0;
	private int cleared = 0;
	
	//each square in the grid - int means empty - other values mean the turn it was placed
	private int[][] field = new int[ROWS][COLS];
	//top row+1 of each column
	//0 means empty
	private int[] top = new int[COLS];
	
	
	//number of next piece
	protected int nextPiece;
	
	
	
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
	
	
	public int[][] getField() {
		return field;
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
	
	// does deep copy of the relevant state information
	public ImprovedState(State s) {
		this(s.getField(), s.getTop(), s.lost, s.getNextPiece(), s.getRowsCleared(), s.getTurnNumber());
	}
	
	//note that this does a deep copy of the parameters. A deep copy is required so that the information can be modified without affecting other objects
	private ImprovedState(int[][] field, int[] top, boolean lost, int nextPiece, int turn, int cleared) {
		for (int row = 0; row < ROWS; row++) {
		System.arraycopy(field[row], 0, this.field[row], 0, COLS);
		}
		
		System.arraycopy(top, 0, this.top, 0, COLS);
		this.nextPiece = nextPiece;
		this.lost = lost;
		this.cleared = cleared;
		this.turn = turn;
	}
	
	//gives legal moves for 
		public int[][] legalMoves() {
			return legalMoves[nextPiece];
		}
		

		//returns a new object representing the state of the game if a particular move was made
		// move should be an index in the legal moves list
				// the returned ImprovedState object ccan then be fed into feature functions
		public ImprovedState tryMove(int move) {
			return tryMove(legalMoves[nextPiece][move]);
		}
		
		//make a move based on an array of orient and slot
		private ImprovedState tryMove(int[] move) {
			return tryMove(move[ORIENT],move[SLOT]);
		}
		
		 
		private ImprovedState tryMove(int orient, int slot) {
			int[][] newField = new int[ROWS][COLS]; 
			
			for (int row = 0; row < ROWS; row++) {
				System.arraycopy(field[row], 0, newField[row], 0, COLS);
				}
				
			int[] newTop = new int[COLS];
				System.arraycopy(top, 0, newTop, 0, COLS);
				int newCleared = cleared;
			//height if the first column makes contact
			int height = newTop[slot]-pBottom[nextPiece][orient][0];
			//for each column beyond the first in the piece
			for(int c = 1; c < pWidth[nextPiece][orient];c++) {
				height = Math.max(height,newTop[slot+c]-pBottom[nextPiece][orient][c]);
			}
			
			//check if game ended
			if(height+pHeight[nextPiece][orient] >= ROWS) {
				return new ImprovedState(newField, newTop, true, nextPiece, turn+1, newCleared);
			}

			
			//for each column in the piece - fill in the appropriate blocks
			for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
				
				//from bottom to top of brick
				for(int h = height+pBottom[nextPiece][orient][i]; h < height+pTop[nextPiece][orient][i]; h++) {
					newField[h][i+slot] = turn;
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
						while(newTop[c]>=1 && newField[top[c]-1][c]==0)	newTop[c]--;
					}
				}
			}

			return new ImprovedState(newField, newTop, false, nextPiece, turn+1, newCleared);
		}

}
