// a convenience function to represent a feature and its associated weight

public class FeatureWeightPair {
	public FeatureFunction feature;
	public float weight;

	FeatureWeightPair(FeatureFunction f, float w) {
		feature = f;
		weight = w;
	}
}
