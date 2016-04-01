/**
 * Created by john on 3/31/2016.
 */
import java.util.ArrayList;
import java.util.concurrent.Callable;

public class Slave implements Callable <Float> {
    private ImprovedState currentState;
    private int move;
    private ArrayList<FeatureWeightPair> features;

    public Slave(ImprovedState currentState, int move, ArrayList<FeatureWeightPair> features) {
        this.currentState = currentState;
        this.features = features;
        this.move = move;
    }

    private float evaluate(ImprovedState s) {
        if(s.hasLost())
            return Float.NEGATIVE_INFINITY;
        float sum = 0.0f;
        for (int i=0; i<features.size(); i++) {
            FeatureWeightPair f = features.get(i);
            sum += f.feature.evaluate(s)*f.weight;
        }
        return sum;
    }

    public Float call() {
        // System.out.println("Slave running: " + move);
        ImprovedState resultingState = currentState.tryMove(move);
        return evaluate(resultingState);
    }
}
