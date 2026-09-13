import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * BiDijkstra re-run driver: re-runs ONLY Bidirectional Dijkstra with the
 * asymmetric (cell-entry) backward relaxation, on exactly the same grids as the
 * main experiment (seed = 100*depth + trial), and emits per-trial rows in the
 * same format as edbt_pertrial.csv. Run AFTER the main run; edbt_pipeline.py merges.
 */
public class BiDRerun {
    public static void main(String[] args) throws Exception {
        int[] DEPTHS = {4, 6, 8, 10, 12};
        int TRIALS   = 30;
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("edbt_bidij_rerun.csv")));
        out.println("depth,trial,alg,N,U,found,cost,deltaS_mean,deltaD,optCost");
        for (int depth : DEPTHS) {
            int gs = 4 * depth;
            for (int trial = 0; trial < TRIALS; trial++) {
                long seed = 100L * depth + trial;
                UnifiedTreeExperimentL.Grid g = new UnifiedTreeExperimentL.Grid(gs, seed);
                UnifiedTreeExperimentL.EDBTMetrics md = UnifiedTreeExperimentL.dijkstraExpand(g);
                double optCost = md.cost;
                UnifiedTreeExperimentL.EDBTMetrics m = UnifiedTreeExperimentL.biDijkstraExpand(g);
                out.println(depth + "," + trial + ",BiDijkstra," + m.nodes + "," + m.useful + ","
                        + (m.found ? 1 : 0)
                        + "," + (!Double.isNaN(m.cost) ? String.format("%.2f", m.cost) : "NA")
                        + "," + String.format("%.6f", m.deltaSMean)
                        + "," + String.format("%.6f", m.deltaD)
                        + "," + (!Double.isNaN(optCost) ? String.format("%.2f", optCost) : "NA"));
                out.flush();
            }
            System.out.println("depth " + depth + " done");
        }
        out.close();
        System.out.println("BiDRerun complete -> edbt_bidij_rerun.csv");
    }
}
