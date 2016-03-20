import java.util.*;
// import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

public class GeneticAlgorithm {
public int populationSize;
ArrayList<ArrayList<FeatureWeightPair>> population;

public GeneticAlgorithm(int populationSize) {
	population = new ArrayList<ArrayList<FeatureWeightPair>>(populationSize);
	
	for (int i = 0; i < populationSize; i++) {
		population.add(generateRandomIndividual());
		/* 
		System.out.print("individual " + i + ":");
		for (FeatureWeightPair f : population.get(i)) {
			System.out.print(f.weight + ", ");
		}
		System.out.print("\n");
		*/
	}
}

private ArrayList<FeatureWeightPair> generateRandomIndividual() {
	ArrayList<FeatureWeightPair> individual = new ArrayList<FeatureWeightPair>();
	// all the feature functions we're using so far contribute negatively to happiness and so should be minimized, 
// 	hence their weights should be negative
	// for example, the presence of holes should decrease happiness
	// not sure if -1.0 is a good lower bound for the initial population
	 	
	individual.add(new FeatureWeightPair(new PlayerSkeleton.AggHeight(), randomFloat(-1.0f, 0.0f)));
	individual.add(new FeatureWeightPair(new PlayerSkeleton.Bumpiness(), randomFloat(-1.0f, 0.0f)));
	individual.add(new FeatureWeightPair(new PlayerSkeleton.MaxHeight(), randomFloat(-1.0f, 0.0f)));
	individual.add(new FeatureWeightPair(new PlayerSkeleton.NumHoles(), randomFloat(-1.0f, 0.0f)));

	return individual;
}

private float randomFloat(float minInclusive, float maxExclusive) {
 return (float) ThreadLocalRandom.current().nextDouble(minInclusive, maxExclusive);
}

	public static void main(String[] args) {
		GeneticAlgorithm ga = new GeneticAlgorithm(100);
		
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



	}

}
