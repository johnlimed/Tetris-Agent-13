// a convenience function to represent a feature and its associated weight

public class FeatureWeightPair {
	public FeatureFunction feature;
	public float weight;
	public boolean increasesHappiness; // used in mutation to determine whether
										// weight values should be positive or
										// negative

	FeatureWeightPair(FeatureFunction f, float w, boolean increasesHappiness) {
		feature = f;
		weight = w;
		this.increasesHappiness = increasesHappiness;
	}
}
