// runs PlayerSkeleton for n random games to obtain information about agent performance for report
// the output to the console should be redirected to a .csv file, and opened in excel for analysis/drawing

import java.util.ArrayList;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class Benchmarker {
	public static void main(String[] args) throws InterruptedException {
		int games = 0;
		Scanner sc = new Scanner(System.in);
		games = sc.nextInt();
		int cores = Runtime.getRuntime().availableProcessors();
		ExecutorService service = Executors.newFixedThreadPool(cores);
		ArrayList<Integer> scores = new ArrayList<>(games);
		ArrayList<Future<Integer>> tasks = new ArrayList<Future<Integer>>(games);

		System.out.println("game number, rows completed");

		for (int i = 0; i < games; i++) {
			PlayerSkeleton p = new PlayerSkeleton();
			ArrayList<FeatureWeightPair> fwPairs = new ArrayList<FeatureWeightPair>();
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.MaxHeight(), -0.87448764f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.MeanHeight(), -0.5054159f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.NumHoles(), -7.5740294f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.RowTransitions(), -1.1404463f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.StdDevHeight(), -1.9600217f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.SumOfPitDepth(), -1.4352853f, false));
			fwPairs.add(new FeatureWeightPair(new PlayerSkeleton.WellSums(), -0.9274553f, false));
			p.setFeatureWeightPairs(fwPairs);
			tasks.add(service.submit(p));
		}

		for (int i = 0; i < games; i++) {
			try {
				int score = tasks.get(i).get();
				// scores.add(score);
				System.out.println(i + ", " + score);
			} catch (Exception e) {
				System.out.println("the future computing a game's score was interupted");
			}
		}

		sc.close();
		service.shutdown();
	}
}
