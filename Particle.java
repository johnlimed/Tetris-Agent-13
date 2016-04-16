
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
