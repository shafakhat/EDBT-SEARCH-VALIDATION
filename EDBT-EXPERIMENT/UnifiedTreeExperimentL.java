import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;
import java.io.*;

// =====================================================================
// EDBT FRAMEWORK - EXPERIMENTAL VALIDATION
// ---------------------------------------------------------------------
// Unified experiment driver for the EDBT article. Twelve search
// algorithms (BFS, Dijkstra, A*, WA*(1.5), WA*(2.0), GBFS, BiDijkstra,
// RBFS, IDA*, Alpha-Beta, NoisyPrune p=0.3 / p=0.6) run on 4d x 4d
// weighted grids (obstacle rate 0.15) at depths d in {4,6,8,10,12},
// 30 trials per configuration (1,800 runs). Per-trial EDBT metrics
// (delta_s, Delta_d, BEI, ABEI, ACP, PCCT, q*) are measured online, a
// per-trial CSV is dumped for independent statistical analysis, and a
// cross-check against the article's printed tables is printed at the end.
// =====================================================================
public class UnifiedTreeExperimentL extends JPanel {

    static final int    NUM_TRIALS    = 30;
    static boolean FAST_MODE = false;
    static final int[]  DEPTHS        = {4, 6, 8, 10, 12};
    static final int    BRANCHING     = 4;     // nominal base b (max grid out-degree; minimax branching)
    static final double OBSTACLE_RATE = 0.15;
    static final int    AB_MAX_DEPTH  = 12;   // never binding at tested depths (d <= 12):
                                              // keeps leaf-parent parity uniform (MIN nodes)
                                              // across all tested depths
    static final double EPSILON       = 0.05;

    static PrintWriter csv = null;

    // ============================================================
    // RESULT RECORD
    // ============================================================
    static class AlgResult {
        String  name;
        int     depth;
        double  meanNodes, stdNodes, meanUseful, meanBEI, solutionRate;
        double  meanDeltaD, meanDeltaS, meanDeltaC, meanABEI, meanQstar;
        double  certDeltaS, acpProduct, gammaMinus, gammaPlus, gammaMu;
        double  meanCostRatio;
        boolean acpSatisfied, pcctSatisfied;

        AlgResult(String n, int d, double mn, double sn,
                  double mu, double mb, double sr,
                  double dd, double ds, double dc, double abei,
                  double qs, double cds, double acp,
                  double gm, double gp, double gmu, double cr) {
            name = n; depth = d;
            meanNodes = mn; stdNodes = sn; meanUseful = mu;
            meanBEI = mb; solutionRate = sr;
            meanDeltaD = dd; meanDeltaS = ds; meanDeltaC = dc;
            meanABEI = abei; meanQstar = qs; certDeltaS = cds;
            acpProduct = acp; gammaMinus = gm; gammaPlus = gp; gammaMu = gmu;
            meanCostRatio = cr;
            acpSatisfied  = Double.isFinite(acp) && Double.isFinite(gm) && (acp >= gm - 1e-9);
            pcctSatisfied = Double.isFinite(ds)  && Double.isFinite(cds) && (ds  >= cds - 1e-9);
        }
    }

    // ============================================================
    // EDBT COMPONENT 1: Complexity Profile
    //   gamma+ = [b^beta/(b^beta-1)] * b^(beta*d)   (Dechter-Pearl-style envelope)
    //   gamma- = b^(beta*d)                          (Knuth-Moore effective-branching floor)
    //   gammaMu= gamma+ * (1-rho)^d                  (density-scaled surrogate)
    //   Coherent parametrisation: beta=1/2 -> alpha-beta (gamma- = b^(d/2) minimal
    //   leaves, gamma+ = 2*b^(d/2) total for b=4); beta=1 -> unpruned tree.
    // ============================================================
    static double clampBeta(double beta) {
        if (!Double.isFinite(beta) || beta < 0.01) return 0.01;
        if (beta > 1.0) return 1.0;
        return beta;
    }

    static double[] complexityProfile(double betaRaw, int d, int b, double rho) {
        double beta   = clampBeta(betaRaw);
        double bBeta  = Math.pow(b, beta);
        double denom  = bBeta - 1.0;
        if (Math.abs(denom) < 1e-10) denom = 1e-10;
        double alpha  = bBeta / denom;
        double bBetaD = Math.pow(b, beta * d);
        double gammaPlus  = alpha * bBetaD;
        double gammaMinus = bBetaD;
        double gammaMu    = alpha * Math.pow(1.0 - rho, d) * bBetaD;
        gammaPlus  = Double.isFinite(gammaPlus)  && gammaPlus  > 0 ? gammaPlus  : 1.0;
        gammaMinus = Double.isFinite(gammaMinus) && gammaMinus > 0 ? gammaMinus : 1.0;
        gammaMu    = Double.isFinite(gammaMu)    && gammaMu    > 0 ? gammaMu    : 1.0;
        return new double[]{gammaPlus, gammaMinus, gammaMu};
    }

    // ============================================================
    // EDBT COMPONENT 2: Accuracy Function
    //   delta_s = min Psi(candidates) / min Psi(retained)  in [0,1]
    //   delta_c = 1 - |retained|/|candidates|              in [0,1)
    // ============================================================
    static double[] accuracyMeasures(double[] psiAll, double[] psiKept) {
        if (psiAll == null || psiAll.length == 0)
            return new double[]{1.0, 0.0};
        if (psiKept == null || psiKept.length == 0)
            return new double[]{-1.0, -1.0};
        double minAll = Double.MAX_VALUE;
        for (double v : psiAll) if (v < minAll) minAll = v;
        double minKept = Double.MAX_VALUE;
        for (double v : psiKept) if (v < minKept) minKept = v;
        // Since Pi(S,.) subseteq S (condition C3), min(kept) >= min(all)
        // always holds, so the accuracy-preserved fraction min(all)/min(kept)
        // is naturally <= 1, needs no clamp, equals 1 exactly when the true
        // best candidate was retained, and is < 1 precisely when it was not.
        double deltaS = (minKept > 1e-12) ? minAll / minKept : 1.0;
        double deltaC = 1.0 - (double) psiKept.length / (double) psiAll.length;
        return new double[]{deltaS, deltaC};
    }

    // Sign-adjusted variant for MAX nodes of minimax trees: the local
    // objective is maximisation, so fidelity is max(all)/max(kept).
    static double accuracyMeasuresMax(double[] psiAll, double[] psiKept) {
        if (psiAll == null || psiAll.length == 0) return 1.0;
        if (psiKept == null || psiKept.length == 0) return -1.0;
        double maxAll = -Double.MAX_VALUE;
        for (double v : psiAll) if (v > maxAll) maxAll = v;
        double maxKept = -Double.MAX_VALUE;
        for (double v : psiKept) if (v > maxKept) maxKept = v;
        return (maxKept > 1e-12) ? maxAll / maxKept : 1.0;
    }

    // ============================================================
    // EDBT COMPONENT 3: Convergence Operator  (Definition 3.4)
    //   gamma-/gamma+ = 1-b^(-beta) is independent of d, so the
    //   tolerance condition is either met at every depth or none:
    //   d*(eps) = 0 if beta >= beta*(eps) = log_b(1/eps), else
    //   unreachable (Double.MAX_VALUE).  q* = DeltaD * (1-eps).
    // ============================================================
    static double[] convergenceOperator(double betaRaw, double deltaD,
                                        double eps, int b) {
        double beta  = clampBeta(betaRaw);
        double dStar;
        if (eps <= 0 || eps >= 1) {
            dStar = Double.MAX_VALUE;
        } else {
            double ratio = 1.0 - Math.pow(b, -beta);   // gamma-/gamma+
            dStar = (ratio >= 1.0 - eps) ? 0.0 : Double.MAX_VALUE;
        }
        double qStar = deltaD * (1.0 - eps);
        return new double[]{dStar, qStar};
    }

    static double betaThreshold(double eps, int b) {
        return Math.log(1.0 / eps) / Math.log(b);
    }

    // ============================================================
    // EDBT: Path-level deltaD (product of per-step sibling ratios
    // along the returned path, recomputed over the FULL neighbour
    // set of every path node)
    // ============================================================
    static double computePathDeltaD(Grid g, int[] par, int start, int goal,
                                    double wt, boolean isAdmissible) {
        if (par == null || (par[goal] == -1 && goal != start)) return 1.0;
        List<Integer> pathNodes = new ArrayList<>();
        int cur = goal;
        Set<Integer> seen = new HashSet<>();
        while (cur != start) {
            if (seen.contains(cur)) return 1.0;
            seen.add(cur);
            pathNodes.add(cur);
            if (par[cur] < 0) return 1.0;
            cur = par[cur];
        }
        pathNodes.add(start);
        Collections.reverse(pathNodes);
        if (pathNodes.size() < 2) return 1.0;
        double deltaD  = 1.0;
        double[] gCost = new double[g.total()];
        Arrays.fill(gCost, Double.MAX_VALUE);
        gCost[start] = 0;
        for (int i = 0; i < pathNodes.size() - 1; i++) {
            int u = pathNodes.get(i);
            int v = pathNodes.get(i + 1);
            int edgeCost = 1;
            for (int[] nb : g.nbrs(u))
                if (g.id(nb[0], nb[1]) == v) { edgeCost = nb[2]; break; }
            gCost[v] = gCost[u] + edgeCost;
        }
        for (int i = 0; i < pathNodes.size() - 1; i++) {
            int u = pathNodes.get(i);
            List<int[]> nbrs = g.nbrs(u);
            if (nbrs.isEmpty()) continue;
            double[] psiAll = new double[nbrs.size()];
            for (int j = 0; j < nbrs.size(); j++) {
                int nid = g.id(nbrs.get(j)[0], nbrs.get(j)[1]);
                psiAll[j] = gCost[u] + nbrs.get(j)[2] + wt * g.h(nid);
            }
            int    nextNode = pathNodes.get(i + 1);
            double gNext    = gCost[nextNode];
            double psiNext  = gNext + wt * g.h(nextNode);
            double[] acc    = accuracyMeasures(psiAll, new double[]{psiNext});
            if (acc[0] >= 0) deltaD *= acc[0];
        }
        return deltaD;
    }

    // Path cost (sum of edge weights) along a parent-pointer path.
    static double pathCost(Grid g, int[] par, int start, int goal) {
        if (par == null || (par[goal] == -1 && goal != start)) return Double.NaN;
        double cost = 0; int cur = goal; Set<Integer> seen = new HashSet<>();
        while (cur != start) {
            if (seen.contains(cur)) return Double.NaN;
            seen.add(cur);
            int p = par[cur];
            if (p < 0) return Double.NaN;
            for (int[] nb : g.nbrs(p))
                if (g.id(nb[0], nb[1]) == cur) { cost += nb[2]; break; }
            cur = p;
        }
        return cost;
    }

    // ============================================================
    // PCCT Certificate (Table-8 scheme of the article)
    //   delta_cert = 1/(1+eps_h):
    //     eps_h = 0     -> admissible evaluation (BFS, Dijkstra, A*,
    //                      BiDijkstra, RBFS, IDA*, NoisyPrune, and
    //                      AlphaBeta whose leaf Psi are exact values)
    //     eps_h = w-1   -> WA*(w)
    //     eps_h = inf   -> GBFS (h-only evaluation; the w->infinity
    //                      limit of the WA* family), certificate 0.
    // ============================================================
    static double pcctEpsH(String algName) {
        if (algName.startsWith("WA*(1.5)")) return 0.5;
        if (algName.startsWith("WA*(2.0)")) return 1.0;
        if (algName.equals("GBFS"))         return Double.POSITIVE_INFINITY;
        return 0.0;
    }

    static double pcctCertificate(String algName) {
        double eh = pcctEpsH(algName);
        if (Double.isInfinite(eh)) return 0.0;
        return 1.0 / (1.0 + eh);
    }

    // ============================================================
    // ABEI and BEI   (log base b = 4)
    // ============================================================
    static double computeABEI(double N, double U, double q, double deltaD, int b) {
        if (N < 1 || U < 1) return 0.0;
        double logNorm = Math.log(N + 1) / Math.log(b);
        if (logNorm < 1e-12) return 0.0;
        return (U * q * deltaD) / (N * logNorm);
    }

    static double computeBEI(double N, double U, double q, int b) {
        if (N < 1 || U < 1) return 0.0;
        double logNorm = Math.log(N + 1) / Math.log(b);
        if (logNorm < 1e-12) return 0.0;
        return (U * q) / (N * logNorm);
    }

    // ============================================================
    // GRID
    // ============================================================
    static class Grid {
        final int R, C;
        final int[][] w;
        final boolean[][] blocked;
        final int sr, sc, gr, gc;
        static final int[] DR = {-1, 1,  0, 0};
        static final int[] DC = { 0, 0, -1, 1};

        Grid(int size, long seed) {
            R = C = size;
            w       = new int[R][C];
            blocked = new boolean[R][C];
            sr = sc = 0; gr = R - 1; gc = C - 1;

            Random rng = new Random(seed);
            for (int r = 0; r < R; r++)
                for (int c = 0; c < C; c++)
                    w[r][c] = 1 + rng.nextInt(10);
            for (int r = 0; r < R; r++)
                for (int c = 0; c < C; c++) {
                    if ((r == sr && c == sc) || (r == gr && c == gc)) continue;
                    if (rng.nextDouble() < OBSTACLE_RATE) blocked[r][c] = true;
                }
            // Guaranteed connectivity: top row and rightmost column are open.
            for (int c = 0; c < C; c++) blocked[0][c]     = false;
            for (int r = 0; r < R; r++) blocked[r][C - 1] = false;
        }

        boolean valid(int r, int c) { return r>=0&&r<R&&c>=0&&c<C&&!blocked[r][c]; }
        int id(int r, int c)  { return r * C + c; }
        int row(int id)       { return id / C; }
        int col(int id)       { return id % C; }
        int total()           { return R * C; }
        int h(int r, int c)   { return Math.abs(r - gr) + Math.abs(c - gc); }
        int h(int id)         { return h(row(id), col(id)); }

        List<int[]> nbrs(int r, int c) {
            List<int[]> out = new ArrayList<>(4);
            for (int d = 0; d < 4; d++) {
                int nr = r + DR[d], nc = c + DC[d];
                if (valid(nr, nc)) out.add(new int[]{nr, nc, w[nr][nc]});
            }
            return out;
        }
        List<int[]> nbrs(int id) { return nbrs(row(id), col(id)); }
    }

    // ============================================================
    // PATH LENGTH
    // ============================================================
    static int pathLengthArr(int[] parent, int start, int goal, int maxLen) {
        if (goal == start) return 1;
        if (parent == null || parent[goal] == -1) return 0;
        int len = 1, cur = goal;
        Set<Integer> vis = new HashSet<>();
        while (cur != start) {
            if (vis.contains(cur) || len > maxLen) return 0;
            vis.add(cur); int p = parent[cur];
            if (p < 0) return 0;
            cur = p; len++;
        }
        return len;
    }

    // ============================================================
    // EDBT METRICS
    // ============================================================
    static class EDBTMetrics {
        long    nodes;
        int     useful;
        double  deltaD;
        double  deltaSMean;
        double  deltaCMean;
        int[]   par;
        boolean found;
        double  cost;      // path cost; NaN if not found / not applicable
    }

    // ============================================================
    // ALGORITHM 1: BFS
    // ============================================================
    static EDBTMetrics bfsExpand(Grid g) {
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        int T = g.total();
        int[]     par = new int[T];     Arrays.fill(par, -1);
        boolean[] vis = new boolean[T];
        Queue<Integer> q = new LinkedList<>();
        q.add(start); vis[start] = true;
        long expanded = 0; boolean found = false;
        while (!q.isEmpty()) {
            int u = q.poll(); expanded++;
            if (u == goal) { found = true; break; }
            for (int[] nb : g.nbrs(u)) {
                int nid = g.id(nb[0], nb[1]);
                if (!vis[nid]) { vis[nid] = true; par[nid] = u; q.add(nid); }
            }
        }
        int useful = found ? pathLengthArr(par, start, goal, T) : 0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = found;
        m.deltaD = 1.0; m.deltaSMean = 1.0; m.deltaCMean = 0.0; m.par = par;
        m.cost = found ? pathCost(g, par, start, goal) : Double.NaN;
        return m;
    }

    // ============================================================
    // ALGORITHM 2: DIJKSTRA
    // ============================================================
    static EDBTMetrics dijkstraExpand(Grid g) {
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        int T = g.total();
        double[] dist = new double[T]; Arrays.fill(dist, Double.MAX_VALUE);
        int[]    par  = new int[T];    Arrays.fill(par, -1);
        dist[start] = 0;
        PriorityQueue<double[]> pq = new PriorityQueue<>(Comparator.comparingDouble(x -> x[0]));
        pq.add(new double[]{0, start});
        boolean[] vis = new boolean[T];
        long expanded = 0; boolean found = false;
        double sumDS = 0, sumDC = 0; int steps = 0;
        while (!pq.isEmpty()) {
            double[] cur = pq.poll(); int u = (int) cur[1];
            if (vis[u]) continue;
            vis[u] = true; expanded++;
            if (u == goal) { found = true; break; }
            List<int[]> nbrs = g.nbrs(u);
            if (!nbrs.isEmpty()) {
                double[] psiAll = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++)
                    psiAll[i] = dist[u] + nbrs.get(i)[2];
                List<Double> kl = new ArrayList<>();
                for (int[] nb : nbrs) {
                    int nid = g.id(nb[0], nb[1]); double nd = dist[u] + nb[2];
                    if (!vis[nid] && nd < dist[nid]) kl.add(nd);
                }
                if (!kl.isEmpty()) {
                    double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                    if (acc[0] >= 0) { sumDS += acc[0]; sumDC += acc[1]; steps++; }
                }
            }
            for (int[] nb : g.nbrs(u)) {
                int nid = g.id(nb[0], nb[1]); double nd = dist[u] + nb[2];
                if (nd < dist[nid]) { dist[nid] = nd; par[nid] = u; pq.add(new double[]{nd, nid}); }
            }
        }
        int useful = found ? pathLengthArr(par, start, goal, T) : 0;
        double deltaD = found ? computePathDeltaD(g, par, start, goal, 1.0, true) : 1.0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = found; m.par = par;
        m.deltaD = deltaD;
        m.deltaSMean = (steps > 0) ? sumDS / steps : 1.0;
        m.deltaCMean = (steps > 0) ? sumDC / steps : 0.0;
        m.cost = found ? pathCost(g, par, start, goal) : Double.NaN;
        return m;
    }

    // ============================================================
    // ALGORITHM 3/4: A* / WA*
    // ============================================================
    static EDBTMetrics astarExpand(Grid g, double wt) {
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        int T = g.total();
        double[] gCost = new double[T]; Arrays.fill(gCost, Double.MAX_VALUE);
        int[]    par   = new int[T];    Arrays.fill(par, -1);
        gCost[start] = 0;
        PriorityQueue<double[]> open = new PriorityQueue<>(Comparator.comparingDouble(x -> x[0]));
        open.add(new double[]{wt * g.h(start), start});
        boolean[] closed = new boolean[T];
        long expanded = 0; boolean found = false;
        double sumDS = 0, sumDC = 0; int steps = 0;
        while (!open.isEmpty()) {
            double[] cur = open.poll(); int u = (int) cur[1];
            if (closed[u]) continue;
            closed[u] = true; expanded++;
            if (u == goal) { found = true; break; }
            List<int[]> nbrs = g.nbrs(u);
            if (!nbrs.isEmpty()) {
                double[] psiAll = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++) {
                    int nid = g.id(nbrs.get(i)[0], nbrs.get(i)[1]);
                    psiAll[i] = gCost[u] + nbrs.get(i)[2] + wt * g.h(nid);
                }
                List<Double> kl = new ArrayList<>();
                for (int[] nb : nbrs) {
                    int nid = g.id(nb[0], nb[1]); double ng = gCost[u] + nb[2];
                    if (!closed[nid] && ng < gCost[nid]) kl.add(ng + wt * g.h(nid));
                }
                if (!kl.isEmpty()) {
                    double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                    if (acc[0] >= 0) { sumDS += acc[0]; sumDC += acc[1]; steps++; }
                }
            }
            for (int[] nb : nbrs) {
                int nid = g.id(nb[0], nb[1]);
                if (closed[nid]) continue;
                double ng = gCost[u] + nb[2];
                if (ng < gCost[nid]) {
                    gCost[nid] = ng; par[nid] = u;
                    open.add(new double[]{ng + wt * g.h(nid), nid});
                }
            }
        }
        int useful = found ? pathLengthArr(par, start, goal, T) : 0;
        double deltaD = found ? computePathDeltaD(g, par, start, goal, wt, wt <= 1.0) : 1.0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = found; m.par = par;
        m.deltaD = deltaD;
        m.deltaSMean = (steps > 0) ? sumDS / steps : 1.0;
        m.deltaCMean = (steps > 0) ? sumDC / steps : 0.0;
        m.cost = found ? pathCost(g, par, start, goal) : Double.NaN;
        return m;
    }

    // ============================================================
    // ALGORITHM 3b: NOISYPRUNE A*  (adversarial (C5) stress test)
    // ------------------------------------------------------------
    // Standard A* (w=1, admissible & consistent Manhattan heuristic)
    // except: at every expansion with >= 2 children, the true
    // Psi-argmin child is excluded from relaxation with fixed
    // probability p -- a direct, dosed violation of (C5) targeting
    // exactly the node whose elimination (C5) forbids.  Because the
    // dropped child is never relaxed, the operator is incomplete:
    // runs may fail to reach the goal, which is itself an
    // informative measurement (violating (C5) sacrifices
    // completeness), so the solution rate is reported.
    // ============================================================
    static EDBTMetrics randomPruneAStarExpand(Grid g, double dropBestProb, long rngSeed) {
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        int T = g.total();
        double[] gCost = new double[T]; Arrays.fill(gCost, Double.MAX_VALUE);
        int[]    par   = new int[T];    Arrays.fill(par, -1);
        gCost[start] = 0;
        PriorityQueue<double[]> open = new PriorityQueue<>(Comparator.comparingDouble(x -> x[0]));
        open.add(new double[]{g.h(start), start});
        boolean[] closed = new boolean[T];
        long expanded = 0; boolean found = false;
        double sumDS = 0, sumDC = 0; int steps = 0;
        Random rng = new Random(rngSeed);
        while (!open.isEmpty()) {
            double[] cur = open.poll(); int u = (int) cur[1];
            if (closed[u]) continue;
            closed[u] = true; expanded++;
            if (u == goal) { found = true; break; }
            List<int[]> nbrs = g.nbrs(u);
            if (nbrs.size() >= 2) {
                double[] psiAll = new double[nbrs.size()];
                int argmin = 0;
                for (int i = 0; i < nbrs.size(); i++) {
                    int nid = g.id(nbrs.get(i)[0], nbrs.get(i)[1]);
                    psiAll[i] = gCost[u] + nbrs.get(i)[2] + g.h(nid);
                    if (psiAll[i] < psiAll[argmin]) argmin = i;
                }
                boolean dropBest = rng.nextDouble() < dropBestProb;
                List<Integer> survive = new ArrayList<>();
                for (int i = 0; i < nbrs.size(); i++)
                    if (!(dropBest && i == argmin)) survive.add(i);
                double[] kl = new double[survive.size()];
                for (int i = 0; i < survive.size(); i++) kl[i] = psiAll[survive.get(i)];
                double[] acc = accuracyMeasures(psiAll, kl);
                if (acc[0] >= 0) { sumDS += acc[0]; sumDC += acc[1]; steps++; }
                for (int i : survive) {
                    int[] nb = nbrs.get(i);
                    int nid = g.id(nb[0], nb[1]);
                    if (closed[nid]) continue;
                    double ng = gCost[u] + nb[2];
                    if (ng < gCost[nid]) {
                        gCost[nid] = ng; par[nid] = u;
                        open.add(new double[]{ng + g.h(nid), nid});
                    }
                }
            } else {
                for (int[] nb : nbrs) {
                    int nid = g.id(nb[0], nb[1]);
                    if (closed[nid]) continue;
                    double ng = gCost[u] + nb[2];
                    if (ng < gCost[nid]) {
                        gCost[nid] = ng; par[nid] = u;
                        open.add(new double[]{ng + g.h(nid), nid});
                    }
                }
            }
        }
        int useful = found ? pathLengthArr(par, start, goal, T) : 0;
        double deltaD = found ? computePathDeltaD(g, par, start, goal, 1.0, true) : 1.0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = found; m.par = par;
        m.deltaD = deltaD;
        m.deltaSMean = (steps > 0) ? sumDS / steps : 1.0;
        m.deltaCMean = (steps > 0) ? sumDC / steps : 0.0;
        m.cost = found ? pathCost(g, par, start, goal) : Double.NaN;
        return m;
    }

    // ============================================================
    // ALGORITHM 5: GBFS
    // ============================================================
    static EDBTMetrics gbfsExpand(Grid g) {
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        int T = g.total();
        int[]     par    = new int[T];     Arrays.fill(par, -1);
        boolean[] vis    = new boolean[T];
        boolean[] inOpen = new boolean[T];
        PriorityQueue<int[]> open = new PriorityQueue<>(Comparator.comparingInt(x -> x[0]));
        open.add(new int[]{g.h(start), start}); inOpen[start] = true;
        long expanded = 0; boolean found = false;
        double sumDS = 0, sumDC = 0; int steps = 0;
        while (!open.isEmpty()) {
            int[] cur = open.poll(); int u = cur[1];
            if (vis[u]) continue;
            vis[u] = true; expanded++;
            if (u == goal) { found = true; break; }
            List<int[]> nbrs = g.nbrs(u);
            if (!nbrs.isEmpty()) {
                double[] psiAll = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++)
                    psiAll[i] = g.h(g.id(nbrs.get(i)[0], nbrs.get(i)[1]));
                List<Double> kl = new ArrayList<>();
                for (int[] nb : nbrs) {
                    int nid = g.id(nb[0], nb[1]);
                    if (!vis[nid]) kl.add((double) g.h(nid));
                }
                if (!kl.isEmpty()) {
                    double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                    if (acc[0] >= 0) { sumDS += acc[0]; sumDC += acc[1]; steps++; }
                }
            }
            for (int[] nb : nbrs) {
                int nid = g.id(nb[0], nb[1]);
                if (!vis[nid] && !inOpen[nid]) {
                    if (par[nid] == -1) par[nid] = u;
                    open.add(new int[]{g.h(nid), nid});
                    inOpen[nid] = true;
                }
            }
        }
        int useful = found ? pathLengthArr(par, start, goal, T) : 0;
        double deltaD = found ? computePathDeltaD(g, par, start, goal, 0.0, false) : 1.0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = found; m.par = par;
        m.deltaD = deltaD;
        m.deltaSMean = (steps > 0) ? sumDS / steps : 1.0;
        m.deltaCMean = (steps > 0) ? sumDC / steps : 0.0;
        m.cost = found ? pathCost(g, par, start, goal) : Double.NaN;
        return m;
    }

    // ============================================================
    // ALGORITHM 6: BIDIRECTIONAL DIJKSTRA
    // ============================================================
    static EDBTMetrics biDijkstraExpand(Grid g) {
        int start = g.id(g.sr,g.sc), goal = g.id(g.gr,g.gc), T = g.total();
        double[] dF = new double[T], dB = new double[T];
        Arrays.fill(dF, Double.MAX_VALUE); Arrays.fill(dB, Double.MAX_VALUE);
        dF[start] = 0; dB[goal] = 0;
        boolean[] vF = new boolean[T], vB = new boolean[T];
        int[] parF = new int[T], parB = new int[T];
        Arrays.fill(parF, -1); Arrays.fill(parB, -1);
        PriorityQueue<double[]> pqF = new PriorityQueue<>(Comparator.comparingDouble(x->x[0]));
        PriorityQueue<double[]> pqB = new PriorityQueue<>(Comparator.comparingDouble(x->x[0]));
        pqF.add(new double[]{0, start}); pqB.add(new double[]{0, goal});
        long expanded = 0; double mu = Double.MAX_VALUE; int meet = -1;
        double sumDS = 0, sumDC = 0; int steps = 0;
        while (true) {
            if (pqF.isEmpty() && pqB.isEmpty()) break;
            double topF = pqF.isEmpty() ? Double.MAX_VALUE : pqF.peek()[0];
            double topB = pqB.isEmpty() ? Double.MAX_VALUE : pqB.peek()[0];
            if (topF + topB >= mu) break;
            if (!pqF.isEmpty() && topF <= topB) {
                double[] cur = pqF.poll(); int u = (int) cur[1];
                if (vF[u]) continue;
                vF[u] = true; expanded++;
                List<int[]> nbrs = g.nbrs(u);
                if (!nbrs.isEmpty()) {
                    double[] psiAll = new double[nbrs.size()];
                    for (int i = 0; i < nbrs.size(); i++) psiAll[i] = dF[u] + nbrs.get(i)[2];
                    List<Double> kl = new ArrayList<>();
                    for (int[] nb : nbrs) {
                        int nid = g.id(nb[0], nb[1]); double nd = dF[u] + nb[2];
                        if (nd < dF[nid]) kl.add(nd);
                    }
                    if (!kl.isEmpty()) {
                        double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                        if (acc[0] >= 0) { sumDS += acc[0]; sumDC += acc[1]; steps++; }
                    }
                }
                for (int[] nb : nbrs) {
                    int nid = g.id(nb[0], nb[1]); double nd = dF[u] + nb[2];
                    if (nd < dF[nid]) { dF[nid] = nd; parF[nid] = u; pqF.add(new double[]{nd, nid}); }
                    // both-tentative meeting check: dF/dB labels are realized walk costs,
                    // so mu stays a valid s-t walk cost under the asymmetric entry-cost convention.
                    if (dB[nid] < Double.MAX_VALUE && dF[nid] + dB[nid] < mu) { mu = dF[nid] + dB[nid]; meet = nid; }
                }
            } else if (!pqB.isEmpty()) {
                double[] cur = pqB.poll(); int u = (int) cur[1];
                if (vB[u]) continue;
                vB[u] = true; expanded++;
                for (int[] nb : g.nbrs(u)) {
                    int nid = g.id(nb[0], nb[1]);
                    // edge costs are cell-ENTRY costs w(cell entered), hence asymmetric.
                    // Walking nid->u enters u, so the true nid->u cost is w(u), not w(nid).
                    int wu = g.w[g.row(u)][g.col(u)];
                    double nd = dB[u] + wu;
                    if (nd < dB[nid]) { dB[nid] = nd; parB[nid] = u; pqB.add(new double[]{nd, nid}); }
                    if (dF[nid] < Double.MAX_VALUE && dF[nid] + dB[nid] < mu) { mu = dF[nid] + dB[nid]; meet = nid; }
                }
            } else break;
        }
        int useful = 0;
        double cost = Double.NaN;
        int[] par = null;
        if (meet >= 0) {
            // build the combined s->t parent chain: forward tree up to meet,
            // then the backward chain meet->...->goal (par[p]=cur when parB[cur]=p).
            par = new int[T];
            for (int i = 0; i < T; i++) par[i] = parF[i];
            int cur = meet; boolean okChain = true;
            while (cur != goal) {
                int p = parB[cur];
                if (p < 0) { okChain = false; break; }
                par[p] = cur; cur = p;
            }
            if (okChain) {
                useful = pathLengthArr(par, start, goal, T);
                cost   = pathCost(g, par, start, goal);   // forward entry-cost convention
            } else { par = null; }
        }
        double deltaD = (par != null && useful > 0)
                ? computePathDeltaD(g, par, start, goal, 1.0, true) : 1.0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = (meet >= 0 && par != null && useful > 0); m.par = par;
        m.deltaD = deltaD;
        m.deltaSMean = (steps > 0) ? sumDS / steps : 1.0;
        m.deltaCMean = (steps > 0) ? sumDC / steps : 0.0;
        m.cost = cost;
        return m;
    }

    // ============================================================
    // ALGORITHM 7: RBFS  (Korf-style recursion with best-g
    // memoisation and an outer f-limit schedule; node cap 2e6)
    // ============================================================
    static long   rbfsCount;
    static double rbfsSumDS;
    static int    rbfsSteps;

    static boolean rbfsFound;
    static double  rbfsBestCost;

    static double rbfsRec(Grid g, int cur, double gCur, double fCur,
                           double fLimit, int goal,
                           Set<Integer> onPath,
                           List<int[]> pathEdges,
                           Map<Integer, Double> bestG) {

        if (rbfsCount >= 2_000_000) {
            rbfsFound = false;
            return Double.MAX_VALUE;
        }
        rbfsCount++;

        if (cur == goal) {
            rbfsFound    = true;
            rbfsBestCost = gCur;
            return 0.0;
        }

        // Candidate set = FULL sibling set with would-be Psi;
        // retained set = generated successors (after cycle guard and
        // best-g pruning), the same convention as IDA*.
        List<int[]> allNbrs = g.nbrs(cur);
        double[] psiCand = new double[allNbrs.size()];
        for (int i = 0; i < allNbrs.size(); i++) {
            int[] nb = allNbrs.get(i);
            int nid  = g.id(nb[0], nb[1]);
            psiCand[i] = gCur + nb[2] + g.h(nid);
        }

        List<double[]> succs = new ArrayList<>();
        for (int[] nb : allNbrs) {
            int    nid      = g.id(nb[0], nb[1]);
            double edgeCost = nb[2];
            double ng       = gCur + edgeCost;

            if (onPath.contains(Integer.valueOf(nid))) continue;

            Double knownG = bestG.get(nid);
            if (knownG != null && ng >= knownG) continue;

            bestG.put(nid, ng);

            double nf = Math.max(ng + g.h(nid), fCur);
            succs.add(new double[]{nid, edgeCost, nf});
        }

        if (succs.isEmpty()) return Double.MAX_VALUE;

        double[] psiKept = succs.stream().mapToDouble(s -> s[2]).toArray();
        double[] acc     = accuracyMeasures(psiCand, psiKept);
        if (acc[0] >= 0) { rbfsSumDS += acc[0]; rbfsSteps++; }

        succs.sort(Comparator.comparingDouble(s -> s[2]));

        int maxIter = succs.size() * 2_000_000;
        int iter    = 0;

        while (iter++ < maxIter) {

            succs.sort(Comparator.comparingDouble(s -> s[2]));
            double[] best = succs.get(0);

            if (best[2] >= Double.MAX_VALUE) return Double.MAX_VALUE;

            if (best[2] > fLimit) return best[2];

            double altF    = succs.size() > 1 ? succs.get(1)[2] : Double.MAX_VALUE;
            int    bestNid = (int) best[0];
            double bestEC  = best[1];

            onPath.add(Integer.valueOf(bestNid));
            pathEdges.add(new int[]{cur, bestNid, (int) bestEC});

            double newLimit = (altF >= Double.MAX_VALUE) ? fLimit : Math.min(fLimit, altF);
            double res = rbfsRec(g, bestNid, gCur + bestEC, best[2],
                                  newLimit, goal, onPath, pathEdges, bestG);

            if (rbfsFound) return 0.0;

            pathEdges.remove(pathEdges.size() - 1);
            onPath.remove(Integer.valueOf(bestNid));

            best[2] = (res >= Double.MAX_VALUE) ? Double.MAX_VALUE : res;
        }

        return Double.MAX_VALUE;
    }

    static EDBTMetrics rbfsExpand(Grid g) {
        rbfsCount = 0; rbfsSumDS = 0.0; rbfsSteps = 0;
        rbfsFound = false; rbfsBestCost = Double.MAX_VALUE;

        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);

        Map<Integer, Double> bestG = new HashMap<>();
        bestG.put(start, 0.0);

        Set<Integer> onPath    = new HashSet<>();
        List<int[]>  pathEdges = new ArrayList<>();
        onPath.add(Integer.valueOf(start));

        double fLimit = g.h(start);
        boolean done  = false;

        for (int iter = 0; iter < 100_000 && !done; iter++) {
            rbfsFound = false;
            pathEdges.clear();
            onPath.clear();
            onPath.add(Integer.valueOf(start));
            bestG.clear();
            bestG.put(start, 0.0);

            double res = rbfsRec(g, start, 0, fLimit, fLimit,
                                  goal, onPath, pathEdges, bestG);

            if (rbfsFound) {
                done = true;
            } else if (res >= Double.MAX_VALUE || rbfsCount >= 2_000_000) {
                break;
            } else {
                fLimit = res;
            }
        }

        boolean found  = rbfsFound;
        int     useful = found ? pathEdges.size() + 1 : 0;

        double deltaD = 1.0;
        double cost   = Double.NaN;
        if (found && !pathEdges.isEmpty()) {
            double gCur = 0;
            for (int[] edge : pathEdges) {
                int u = edge[0], v = edge[1], ec = edge[2];
                List<int[]> nbrs = g.nbrs(u);
                double[] psiAll  = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++) {
                    int nid = g.id(nbrs.get(i)[0], nbrs.get(i)[1]);
                    psiAll[i] = gCur + nbrs.get(i)[2] + g.h(nid);
                }
                double gNext    = gCur + ec;
                double[] psiKpt = new double[]{gNext + g.h(v)};
                double[] a      = accuracyMeasures(psiAll, psiKpt);
                if (a[0] >= 0) deltaD *= a[0];
                gCur = gNext;
            }
            cost = gCur;
        }

        EDBTMetrics m = new EDBTMetrics();
        m.nodes      = rbfsCount;
        m.useful     = useful;
        m.found      = found;
        m.par        = null;
        m.deltaD     = deltaD;
        m.deltaSMean = (rbfsSteps > 0) ? rbfsSumDS / rbfsSteps : 1.0;
        m.deltaCMean = 0.0;
        m.cost       = cost;
        return m;
    }

    // ============================================================
    // ALGORITHM 8: IDA*  (iterative deepening, on-path cycle guard,
    // node cap 1e7)
    // ============================================================
    static long   idaCount;
    static double idaSumDS;
    static int    idaStepsCount;

    static double idaSearch(Grid g, int cur, double gCur, double thresh,
                             int goal, Set<Integer> onPath, List<int[]> pathEdges) {
        double f = gCur + g.h(cur);
        if (f > thresh) return f;
        idaCount++;
        if (idaCount > 10_000_000) return Double.MAX_VALUE;
        if (cur == goal) return -gCur;

        List<int[]> nbrs = g.nbrs(cur);
        nbrs.sort(Comparator.comparingInt(nb -> (int)(gCur + nb[2]) + g.h(nb[0], nb[1])));

        if (!nbrs.isEmpty()) {
            double[] psiAll = new double[nbrs.size()];
            for (int i = 0; i < nbrs.size(); i++) {
                int nid = g.id(nbrs.get(i)[0], nbrs.get(i)[1]);
                psiAll[i] = gCur + nbrs.get(i)[2] + g.h(nid);
            }
            List<Double> kl = new ArrayList<>();
            for (int[] nb : nbrs) {
                int nid = g.id(nb[0], nb[1]);
                if (!onPath.contains(Integer.valueOf(nid)))
                    kl.add(gCur + nb[2] + (double) g.h(nid));
            }
            if (!kl.isEmpty()) {
                double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                if (acc[0] >= 0) { idaSumDS += acc[0]; idaStepsCount++; }
            }
        }

        double minExc = Double.MAX_VALUE;
        for (int[] nb : nbrs) {
            int nid = g.id(nb[0], nb[1]);
            if (onPath.contains(Integer.valueOf(nid))) continue;
            onPath.add(Integer.valueOf(nid));
            pathEdges.add(new int[]{cur, nid, nb[2]});
            double res = idaSearch(g, nid, gCur + nb[2], thresh, goal, onPath, pathEdges);
            if (res < 0) return res;
            pathEdges.remove(pathEdges.size() - 1);
            onPath.remove(Integer.valueOf(nid));
            if (res < minExc) minExc = res;
        }
        return minExc;
    }

    static EDBTMetrics idaExpand(Grid g) {
        idaCount = 0; idaSumDS = 0.0; idaStepsCount = 0;
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        double thresh = g.h(start);
        Set<Integer> onPath    = new HashSet<>();
        List<int[]>  pathEdges = new ArrayList<>();
        onPath.add(Integer.valueOf(start));
        boolean found = false;
        for (int iter = 0; iter < 100_000; iter++) {
            double res = idaSearch(g, start, 0, thresh, goal, onPath, pathEdges);
            if (res < 0)                   { found = true; break; }
            if (res == Double.MAX_VALUE)   break;
            thresh = res;
        }
        int useful = found ? pathEdges.size() + 1 : 0;
        double deltaD = 1.0;
        double cost = Double.NaN;
        if (found && !pathEdges.isEmpty()) {
            double gCur = 0;
            for (int[] edge : pathEdges) {
                int u = edge[0], v = edge[1], ec = edge[2];
                List<int[]> nbrs = g.nbrs(u);
                double[] psiAll  = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++)
                    psiAll[i] = gCur + nbrs.get(i)[2]
                              + g.h(g.id(nbrs.get(i)[0], nbrs.get(i)[1]));
                double gNext = gCur + ec;
                double[] a   = accuracyMeasures(psiAll, new double[]{gNext + g.h(v)});
                if (a[0] >= 0) deltaD *= a[0];
                gCur = gNext;
            }
            cost = gCur;
        }
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = idaCount; m.useful = useful; m.found = found; m.par = null;
        m.deltaD = deltaD;
        m.deltaSMean = (idaStepsCount > 0) ? idaSumDS / idaStepsCount : 1.0;
        m.deltaCMean = 0.0;
        m.cost = cost;
        return m;
    }

    // ============================================================
    // ALGORITHM 9: ALPHA-BETA  (b=4 uniform minimax tree, leaf values
    // uniform in {0..100}; move ordering is best-first at MAX nodes
    // and worst-first (reverse-sorted) at MIN nodes, guaranteeing
    // genuinely imperfect global ordering; effective depth capped at
    // AB_MAX_DEPTH = 9).
    // ============================================================
    static long   abCount;
    static double abSumDS;
    static int    abSteps;
    static double abSumDMax;     // separate accumulators for MAX-level fidelity
    static int    abStepsMax;

    static int alphaBeta(int[] tree, int node, int depth,
                          int alpha, int beta, boolean maxNode,
                          int b, int maxDepth, int firstLeaf) {
        abCount++;
        if (node >= firstLeaf || depth >= maxDepth)
            return tree[Math.min(node, tree.length - 1)];
        int[] children = childrenOf(node, b, tree.length);

        // Online fidelity is measured only where Psi is
        // materialised online: at nodes whose children are leaves.
        boolean leafParent = (children[0] >= firstLeaf);
        double[] psiAll = null;
        if (leafParent) {
            psiAll = new double[children.length];
            for (int i = 0; i < children.length; i++)
                psiAll[i] = tree[Math.min(children[i], tree.length - 1)];
        }
        List<Double> keptPsi = new ArrayList<>();

        if (maxNode) {
            Integer[] box = new Integer[children.length];
            for (int i = 0; i < children.length; i++) box[i] = children[i];
            Arrays.sort(box, (a2, b2) -> Integer.compare(
                    tree[Math.min(b2, tree.length-1)], tree[Math.min(a2, tree.length-1)]));
            for (int i = 0; i < children.length; i++) children[i] = box[i];
            int val = Integer.MIN_VALUE;
            for (int ch : children) {
                int cv = alphaBeta(tree, ch, depth+1, alpha, beta, false, b, maxDepth, firstLeaf);
                keptPsi.add((double) cv);
                val = Math.max(val, cv); alpha = Math.max(alpha, val);
                if (beta <= alpha) break;
            }
            if (leafParent) {
                // sign-adjusted fidelity at MAX nodes: max(cand)/max(kept)
                double dmax = accuracyMeasuresMax(psiAll, keptPsi.stream().mapToDouble(x->x).toArray());
                if (dmax >= 0) { abSumDMax += dmax; abStepsMax++; }
            }
            return val;
        } else {
            Integer[] box = new Integer[children.length];
            for (int i = 0; i < children.length; i++) box[i] = children[i];
            Arrays.sort(box, (a2, b2) -> Integer.compare(
                    tree[Math.min(b2, tree.length-1)], tree[Math.min(a2, tree.length-1)]));
            for (int i = 0; i < children.length; i++) children[i] = box[i];
            int val = Integer.MAX_VALUE;
            for (int ch : children) {
                int cv = alphaBeta(tree, ch, depth+1, alpha, beta, true, b, maxDepth, firstLeaf);
                keptPsi.add((double) cv);
                val = Math.min(val, cv);
                beta = Math.min(beta, val);
                if (beta <= alpha) break;
            }
            if (leafParent) {
                double[] acc = accuracyMeasures(psiAll, keptPsi.stream().mapToDouble(x->x).toArray());
                if (acc[0] >= 0) { abSumDS += acc[0]; abSteps++; }
            }
            return val;
        }
    }

    static int[] childrenOf(int node, int b, int treeLen) {
        int[] ch = new int[b];
        for (int i = 0; i < b; i++) {
            ch[i] = b * node + 1 + i;
            if (ch[i] >= treeLen) ch[i] = treeLen - 1;
        }
        return ch;
    }

    // Note: no path-level Delta_d is computed for AlphaBeta.
    // The returned principal variation of an exact minimax search is
    // locally argoptimal BY CONSTRUCTION at every level (minimax
    // backing up exact leaf values), so any exact-value path
    // measurement is trivially 1; the informative fidelity signal
    // for adversarial search lives at the node level (delta_s), where
    // cutoffs under imperfect ordering measurably discard the
    // Psi-optimal sibling.  Delta_d is therefore reported as 1.0 by
    // convention, exactly as for BFS, and this convention is stated
    // explicitly in the article.

    static EDBTMetrics abExpand(int b, int depth, long seed) {
        abCount = 0; abSumDS = 0.0; abSteps = 0; abSumDMax = 0.0; abStepsMax = 0;
        Random rng = new Random(seed);
        int safeDepth = Math.min(depth, AB_MAX_DEPTH);
        long size = 1, acc = 1;
        for (int i = 0; i < safeDepth; i++) {
            acc *= b;
            if (acc > 20_000_000L) { safeDepth = i; break; }
            size += acc;
        }
        int[] tree      = new int[(int) size];
        int   firstLeaf = (int) (size - acc);
        for (int i = firstLeaf; i < tree.length; i++) tree[i] = rng.nextInt(101);
        alphaBeta(tree, 0, 0, Integer.MIN_VALUE, Integer.MAX_VALUE,
                  true, b, safeDepth, firstLeaf);

        // delta_s = mean over leaf-parent measurements (MIN nodes via
        // min-ratio, MAX nodes via max-ratio).
        double deltaS;
        if (abSteps + abStepsMax > 0)
            deltaS = (abSumDS + abSumDMax) / (abSteps + abStepsMax);
        else deltaS = 1.0;
        double deltaSMin = (abSteps  > 0) ? abSumDS  / abSteps  : 1.0;
        double deltaSMax = (abStepsMax > 0) ? abSumDMax / abStepsMax : 1.0;
        System.out.printf("   [AlphaBeta detail] ds(MIN leaf-parents)=%.3f (n=%d), ds(MAX leaf-parents)=%.3f (n=%d)%n",
                deltaSMin, abSteps, deltaSMax, abStepsMax);

        double deltaD = 1.0;   // by construction; see note above

        EDBTMetrics m = new EDBTMetrics();
        m.nodes = abCount; m.useful = safeDepth + 1; m.found = true; m.par = null;
        m.deltaD = deltaD;
        m.deltaSMean = deltaS;
        m.deltaCMean = 0.0;
        m.cost = Double.NaN;   // minimax value; not comparable with grid path costs
        return m;
    }

    // ============================================================
    // BETA ESTIMATION  (Eq. (4): incremental least-squares slope of
    // log N on d*log b, cold-start default beta = 0.5)
    // ============================================================
    static double[] fitBeta(List<Integer> depths, List<Double> meanNodes, int b) {
        int n = depths.size();
        if (n < 2) return new double[]{0.5, 1.0, 0.0};
        double logb = Math.log(b);
        double[] X = new double[n], Y = new double[n];
        for (int i = 0; i < n; i++) {
            X[i] = depths.get(i) * logb;
            Y[i] = Math.log(Math.max(meanNodes.get(i), 1.0));
        }
        double Xm = 0, Ym = 0;
        for (int i = 0; i < n; i++) { Xm += X[i]; Ym += Y[i]; }
        Xm /= n; Ym /= n;
        double sXX = 0, sXY = 0, sYY = 0;
        for (int i = 0; i < n; i++) {
            sXX += (X[i]-Xm)*(X[i]-Xm);
            sXY += (X[i]-Xm)*(Y[i]-Ym);
            sYY += (Y[i]-Ym)*(Y[i]-Ym);
        }
        double beta  = (sXX < 1e-12) ? 0.5 : sXY / sXX;
        beta         = clampBeta(beta);
        double a     = Ym - beta * Xm;
        double alpha = Math.exp(a);
        double r2    = (sYY < 1e-12) ? 1.0 : (sXY * sXY) / (sXX * sYY);
        return new double[]{beta, alpha, r2};
    }

    // ============================================================
    // EXPERIMENT RUNNER
    // ============================================================
    static List<AlgResult>           results   = new ArrayList<>();
    static Map<String,List<Double>>  algMeans  = new LinkedHashMap<>();
    static Map<String,List<Integer>> algDepths = new LinkedHashMap<>();
    static Map<String,List<Double>>  vizNodes  = new LinkedHashMap<>();

    static final String[] ALL_ALGS = {
        "BFS","Dijkstra","A*","WA*(1.5)","WA*(2.0)","GBFS",
        "BiDijkstra","RBFS","IDA*","AlphaBeta",
        "NoisyPrune(p=0.3)","NoisyPrune(p=0.6)"
    };

    static void runExperiments() {
        System.out.println("=".repeat(80));
        System.out.println(" EDBT FRAMEWORK - EXPERIMENTAL VALIDATION");
        System.out.printf(" Grid=4d*4d | rho=%.2f | eps=%.2f | Trials=%d | b=%d%n%n",
                OBSTACLE_RATE, EPSILON, NUM_TRIALS, BRANCHING);

        for (String a : ALL_ALGS) {
            algMeans.put(a, new ArrayList<>());
            algDepths.put(a, new ArrayList<>());
            vizNodes.put(a, new ArrayList<>());
        }

        for (int depth : DEPTHS) {
            int gs = 4 * depth;
            Map<String,Long>   sN  = new LinkedHashMap<>(), sN2 = new LinkedHashMap<>();
            Map<String,Long>   sU  = new LinkedHashMap<>(), sF  = new LinkedHashMap<>();
            Map<String,Double> sDD = new LinkedHashMap<>(), sDS = new LinkedHashMap<>(),
                               sDC = new LinkedHashMap<>(), sDDs = new LinkedHashMap<>(),
                               sCR = new LinkedHashMap<>();
            for (String a : ALL_ALGS) {
                sN.put(a,0L); sN2.put(a,0L); sU.put(a,0L); sF.put(a,0L);
                sDD.put(a,0.0); sDS.put(a,0.0); sDC.put(a,0.0);
                sDDs.put(a,0.0); sCR.put(a,0.0);
            }

            for (int trial = 0; trial < NUM_TRIALS; trial++) {
                long seed = 100L * depth + trial;
                Grid g = new Grid(gs, seed);
                Map<String, EDBTMetrics> tr = new LinkedHashMap<>();
                tr.put("BFS",           bfsExpand(g));
                tr.put("Dijkstra",      dijkstraExpand(g));
                // optimal reference cost from Dijkstra
                double optCost = tr.get("Dijkstra").cost;
                tr.put("A*",             astarExpand(g, 1.0));
                tr.put("WA*(1.5)",       astarExpand(g, 1.5));
                tr.put("WA*(2.0)",       astarExpand(g, 2.0));
                tr.put("GBFS",           gbfsExpand(g));
                tr.put("BiDijkstra",     biDijkstraExpand(g));
                if (!FAST_MODE) {
                    tr.put("RBFS",          rbfsExpand(g));
                    tr.put("IDA*",          idaExpand(g));
                }
                tr.put("AlphaBeta",      abExpand(BRANCHING, depth, seed));
                tr.put("NoisyPrune(p=0.3)", randomPruneAStarExpand(g, 0.3, 900_000_000L + seed));
                tr.put("NoisyPrune(p=0.6)", randomPruneAStarExpand(g, 0.6, 800_000_000L + seed));
                for (String alg : tr.keySet()) {
                    EDBTMetrics m = tr.get(alg);
                    acc(sN,sN2,sU,sF,sDD,sDS,sDC,sDDs,sCR,alg,m);
                    // per-trial realised-suboptimality ratio C/C* (successful runs)
                    if (m.found && !Double.isNaN(m.cost) && m.cost > 0
                            && !alg.equals("AlphaBeta") && !Double.isNaN(optCost) && optCost > 0)
                        sCR.merge(alg, m.cost / optCost, Double::sum);
                    pendingRows.add(alg + "," + m.nodes + "," + m.useful + "," + (m.found?1:0)
                            + "," + (!Double.isNaN(m.cost) ? String.format("%.2f", m.cost) : "NA")
                            + "," + String.format("%.6f", m.deltaSMean)
                            + "," + String.format("%.6f", m.deltaD)
                            + "," + (!Double.isNaN(optCost) ? String.format("%.2f", optCost) : "NA"));
                }
                dumpTrial(depth, trial);
            }

            System.out.printf("Depth=%2d (grid=%dx%d):%n", depth, gs, gs);
            for (String alg : ALL_ALGS) {
                double mn  = sN.get(alg)  / (double) NUM_TRIALS;
                double mu  = sU.get(alg)  / (double) NUM_TRIALS;
                double sr  = sF.get(alg)  / (double) NUM_TRIALS;
                double var = sN2.get(alg) / (double) NUM_TRIALS - mn * mn;
                double std = Math.sqrt(Math.max(var, 0));
                // Delta_d averaged over successful runs only
                double dd  = (sF.get(alg) > 0) ? sDDs.get(alg) / sF.get(alg) : 1.0;
                double ds  = sDS.get(alg) / (double) NUM_TRIALS;
                double dc  = sDC.get(alg) / (double) NUM_TRIALS;
                double cr  = (sF.get(alg) > 0 && sCR.get(alg) > 0 && !alg.equals("AlphaBeta"))
                           ? sCR.get(alg) / sF.get(alg) : Double.NaN;
                double q   = alg.startsWith("WA*(1.5)") ? 1.0/1.5
                           : alg.startsWith("WA*(2.0)") ? 1.0/2.0 : 1.0;
                double bei  = computeBEI(mn, mu, q, BRANCHING);
                double abei = computeABEI(mn, mu, q, dd, BRANCHING);

                List<Double>  mns = algMeans.get(alg);
                List<Integer> ds2 = algDepths.get(alg);
                mns.add(mn); ds2.add(depth);

                double[] fit  = fitBeta(ds2, mns, BRANCHING);
                double   beta = fit[0];
                double[] gam  = complexityProfile(beta, depth, BRANCHING, OBSTACLE_RATE);
                double gamP = gam[0], gamM = gam[1], gamMu = gam[2];
                double acp  = mn * dd;
                double cert = pcctCertificate(alg);
                double[] om = convergenceOperator(beta, dd, EPSILON, BRANCHING);

                AlgResult ar = new AlgResult(alg, depth, mn, std, mu, bei, sr,
                        dd, ds, dc, abei, om[1], cert, acp, gamM, gamP, gamMu, cr);
                results.add(ar);
                vizNodes.get(alg).add(mn);

                System.out.printf(
                    "  %-16s N=%9.1f+/-%7.1f SR=%4.1f%% | Dd=%.3f | ds=%.3f | " +
                    "BEI=%.4f | ABEI=%.4f | ACP=%.1f vs g-=%.1f[%s] | q*=%.3f | C/C*=%.3f%n",
                    alg, mn, std, sr*100.0, dd, ds, bei, abei, acp, gamM,
                    ar.acpSatisfied ? "OK" : "VIOLATED", om[1], cr);
            }
            System.out.println();
            System.out.flush();
        }

        printBetaTable();
        printComplexityProfileTable();
        printACPTable();
        printPCCTTable();
        printABEIvsBEI();
        printConvergenceTable();
        printBEITable();
        printHeldOut();
        printRanking();
    }

    // per-trial records for the CSV dump
    static java.util.List<String> trialRows = new ArrayList<>();

    static void dumpTrial(int depth, int trial) {
        for (String row : pendingRows) {
            trialRows.add(depth + "," + trial + "," + row);
        }
        pendingRows.clear();
    }
    static java.util.List<String> pendingRows = new ArrayList<>();

    static void acc(Map<String,Long> sN, Map<String,Long> sN2,
                    Map<String,Long> sU, Map<String,Long> sF,
                    Map<String,Double> sDD, Map<String,Double> sDS,
                    Map<String,Double> sDC, Map<String,Double> sDDs,
                    Map<String,Double> sCR,
                    String alg, EDBTMetrics m) {
        sN.merge(alg, m.nodes, Long::sum);
        sN2.merge(alg, m.nodes * m.nodes, Long::sum);
        sU.merge(alg, (long) m.useful, Long::sum);
        if (m.found) sF.merge(alg, 1L, Long::sum);
        sDD.merge(alg, m.deltaD, Double::sum);
        sDS.merge(alg, m.deltaSMean, Double::sum);
        sDC.merge(alg, m.deltaCMean, Double::sum);
        if (m.found) sDDs.merge(alg, m.deltaD, Double::sum);
    }

    // ============================================================
    // TABLES
    // ============================================================
    static void printBetaTable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: ESTIMATED BETA PARAMETERS (b=" + BRANCHING + ")");
        System.out.printf(" %-16s | %8s | %10s | %6s%n","Algorithm","beta","alpha","R2");
        System.out.println(" "+"-".repeat(48));
        for (String alg : ALL_ALGS) {
            List<Double> mn = algMeans.get(alg); List<Integer> d = algDepths.get(alg);
            if (mn.isEmpty()) continue;
            double[] fit = fitBeta(d, mn, BRANCHING);
            System.out.printf(" %-16s | %8.4f | %10.4f | %6.4f%n", alg, fit[0], fit[1], fit[2]);
        }
        System.out.println();
    }

    static void printComplexityProfileTable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: COMPLEXITY PROFILE BOUNDS");
        for (String alg : ALL_ALGS) {
            List<Double> mn = algMeans.get(alg); List<Integer> d = algDepths.get(alg);
            if (mn.isEmpty()) continue;
            double beta = fitBeta(d, mn, BRANCHING)[0];
            System.out.printf(" %s (beta=%.4f):%n", alg, beta);
            System.out.printf("  %4s|%10s|%10s|%10s|%10s|%s%n","d","g-","gMu","Nbar","g+","In[g-,g+]?");
            boolean allOk = true;
            for (AlgResult r : results) {
                if (!r.name.equals(alg)) continue;
                boolean ok = r.meanNodes >= r.gammaMinus - 1e-3
                          && r.meanNodes <= r.gammaPlus + r.gammaPlus * 0.05;
                if (!ok) allOk = false;
                System.out.printf("  %4d|%10.1f|%10.1f|%10.1f|%10.1f|%s%n",
                        r.depth, r.gammaMinus, r.gammaMu, r.meanNodes, r.gammaPlus, ok?"OK":"no");
            }
            System.out.printf("  -> All OK: %s%n%n", allOk ? "YES" : "NO");
        }
    }

    static void printACPTable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: ACP BOUND  Nbar*Dd >= g-(beta,d)   [ratios]");
        System.out.printf(" %-16s|%4s|%10s|%8s|%10s|%10s|%8s|%s%n",
                "Algorithm","d","Nbar","Dd","Nbar*Dd","g-","ratio","OK?");
        System.out.println(" "+"-".repeat(80));
        for (AlgResult r : results)
            System.out.printf(" %-16s|%4d|%10.1f|%8.3f|%10.2f|%10.2f|%8.2f|%s%n",
                    r.name, r.depth, r.meanNodes, r.meanDeltaD,
                    r.acpProduct, r.gammaMinus,
                    r.acpProduct / r.gammaMinus, r.acpSatisfied?"OK":"VIOLATED");
        System.out.println();
    }

    static void printPCCTTable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: PCCT  ds_obs >= ds_cert = 1/(1+eps_h)");
        System.out.printf(" %-16s|%4s|%9s|%9s|%9s|%9s|%s%n",
                "Algorithm","d","eps_h","ds_cert","ds_obs","Dd","OK?");
        System.out.println(" "+"-".repeat(70));
        for (AlgResult r : results) {
            double eh = pcctEpsH(r.name);
            System.out.printf(" %-16s|%4d|%9s|%9.3f|%9.3f|%9.3f|%s%n",
                    r.name, r.depth, Double.isInfinite(eh)?"inf":String.format("%.2f",eh),
                    r.certDeltaS, r.meanDeltaS, r.meanDeltaD, r.pcctSatisfied?"OK":"FAIL");
        }
        System.out.println();
    }

    static void printABEIvsBEI() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: ABEI vs BEI at d=10");
        System.out.printf(" %-16s|%9s|%7s|%8s|%8s|%8s|%s%n",
                "Algorithm","Nbar","Dd","BEI","ABEI","ratio","ABEI vs BEI");
        System.out.println(" "+"-".repeat(70));
        results.stream().filter(r -> r.depth == 10)
               .sorted(Comparator.comparingDouble((AlgResult r) -> r.meanABEI).reversed())
               .forEach(r -> System.out.printf(
                   " %-16s|%9.1f|%7.3f|%8.4f|%8.4f|%8.1f|%s%n",
                   r.name, r.meanNodes, r.meanDeltaD, r.meanBEI, r.meanABEI,
                   r.meanBEI > 0 ? r.meanBEI / r.meanABEI : 0.0,
                   r.meanABEI > r.meanBEI ? "UP" : r.meanABEI < r.meanBEI ? "DOWN" : "="));
        System.out.println();
    }

    static void printConvergenceTable() {
        System.out.println("=".repeat(80));
        System.out.printf(" TABLE: CONVERGENCE OPERATOR Omega(eps=%.2f) at d=10   beta*(eps)=log_b(1/eps)=%.3f%n",
                EPSILON, betaThreshold(EPSILON, BRANCHING));
        System.out.printf(" %-16s|%8s|%7s|%10s|%8s|%s%n","Algorithm","beta","Dd","beta>=b*?","d*","q*");
        System.out.println(" "+"-".repeat(64));
        for (AlgResult r : results) {
            if (r.depth != 10) continue;
            List<Double> mn = algMeans.get(r.name); List<Integer> d = algDepths.get(r.name);
            double beta = fitBeta(d, mn, BRANCHING)[0];
            double[] om = convergenceOperator(beta, r.meanDeltaD, EPSILON, BRANCHING);
            System.out.printf(" %-16s|%8.4f|%7.3f|%10s|%8s|%8.4f%n",
                    r.name, beta, r.meanDeltaD,
                    beta >= betaThreshold(EPSILON, BRANCHING) ? "yes" : "no",
                    om[0] == Double.MAX_VALUE ? "inf" : String.format("%.0f", om[0]),
                    om[1]);
        }
        System.out.println();
    }

    static void printBEITable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: BEI AND ABEI ACROSS DEPTHS");
        System.out.printf(" %-16s","Algorithm");
        for (int d : DEPTHS) System.out.printf("|d=%-2d  BEI    ABEI  ",d);
        System.out.println();
        System.out.println(" "+"-".repeat(16 + DEPTHS.length * 20));
        for (String alg : ALL_ALGS) {
            System.out.printf(" %-16s", alg);
            for (AlgResult r : results)
                if (r.name.equals(alg))
                    System.out.printf("|%7.4f %7.4f ", r.meanBEI, r.meanABEI);
            System.out.println();
        }
        System.out.println();
    }

    static void printHeldOut() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: HELD-OUT PREDICTION ERROR (train d=4,6,8; test d=10,12)");
        for (String alg : ALL_ALGS) {
            List<Double> mn = algMeans.get(alg); List<Integer> d = algDepths.get(alg);
            if (mn.size() < 5) continue;
            double[] fit     = fitBeta(d.subList(0,3), mn.subList(0,3), BRANCHING);
            double totalErr  = 0; int cnt = 0;
            for (int i = 3; i < d.size(); i++) {
                double pred = fit[1] * Math.pow(BRANCHING, fit[0] * d.get(i));
                double act  = mn.get(i);
                if (act > 0) { totalErr += Math.abs(pred - act) / act; cnt++; }
            }
            System.out.printf(" %-16s | MAPE = %7.2f%%%n",
                    alg, cnt > 0 ? totalErr / cnt * 100 : 0);
        }
        System.out.println();
    }

    static void printRanking() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: RANKING AT d=12 (by ABEI)");
        results.stream().filter(r -> r.depth == 12)
               .sorted(Comparator.comparingDouble((AlgResult r) -> r.meanABEI).reversed())
               .forEach(new java.util.function.Consumer<AlgResult>() {
                   int rank = 1;
                   public void accept(AlgResult r) {
                       System.out.printf(
                           " #%2d %-16s N=%9.1f | SR=%4.1f%% | BEI=%.4f | ABEI=%.4f | Dd=%.3f | q*=%.3f | C/C*=%.3f%n",
                           rank++, r.name, r.meanNodes, r.solutionRate*100.0, r.meanBEI,
                           r.meanABEI, r.meanDeltaD, r.meanQstar,
                           Double.isNaN(r.meanCostRatio) ? -1 : r.meanCostRatio);
                   }
               });
        System.out.println();
    }

    // ============================================================
    // ARTICLE CROSS-CHECK (validation only; does not affect computation)
    // ============================================================
    static final double[][] ARTICLE_N = {
        {222.3,220.2,216.9,212.3,203.1,37.0,155.6,47617.4,8636986.9,137.2,215.5,182.4},
        {497.6,495.6,491.3,483.8,466.7,55.8,359.0,310305.7,10000020.9,1082.5,488.0,413.2},
        {879.6,878.0,874.7,866.4,846.1,79.0,639.3,1263091.2,10000020.2,8681.2,868.9,735.7},
        {1366.7,1364.8,1360.8,1353.9,1325.2,95.1,1008.5,2000000.0,10000020.5,22294.1,1351.8,1144.6},
        {1970.3,1968.3,1965.4,1958.6,1930.0,114.7,1464.6,2000000.0,10000020.9,22218.8,1952.4,1653.3}
    };
    static final double[][] ARTICLE_DD = {
        {1.000,.826,.829,.889,.922,.351,.884,1.000,1.000,.692,.799,.902},
        {1.000,.751,.754,.839,.885,.208,.832,1.000,1.000,.575,.715,.857},
        {1.000,.682,.682,.791,.850,.123,.782,1.000,1.000,.478,.639,.814},
        {1.000,.620,.618,.746,.816,.073,.735,1.000,1.000,.436,.571,.773},
        {1.000,.564,.561,.704,.784,.043,.692,1.000,1.000,.436,.511,.734}
    };
    static final String[] ARTICLE_ALGS = {
        "BFS","Dijkstra","A*","WA*(1.5)","WA*(2.0)","GBFS",
        "BiDijkstra","RBFS","IDA*","AlphaBeta","NoisyPrune(p=0.3)","NoisyPrune(p=0.6)"
    };

    static void printArticleCrossCheck() {
        System.out.println("============================================================");
        System.out.println(" ARTICLE CROSS-CHECK vs the article's printed tables");
        System.out.println(" (N: expected to PASS;  DeltaD: expected to differ where the");
                System.out.println("  printed cell was rounded or re-measured)");
        System.out.println("============================================================");
        int n=0, np=0, dpass=0, dcount=0;
        for (AlgResult r : results) {
            int ai = -1;
            for (int j=0;j<ARTICLE_ALGS.length;j++) if (ARTICLE_ALGS[j].equals(r.name)) { ai=j; break; }
            int di=-1;
            for (int j=0;j<DEPTHS.length;j++) if (DEPTHS[j]==r.depth) { di=j; break; }
            if (ai<0 || di<0) continue;
            double nerr=Math.abs(r.meanNodes-ARTICLE_N[di][ai]);
            double derr=Math.abs(r.meanDeltaD-ARTICLE_DD[di][ai]);
            boolean npass=nerr<=0.15;
            boolean dpassOne=derr<=0.005;
            n++; dcount++;
            if(npass)np++; if(dpassOne)dpass++;
            System.out.printf("%-18s d=%2d | N %9.2f vs %9.2f [%s] | Dd %.3f vs %.3f [%s]%n",
                r.name,r.depth,r.meanNodes,ARTICLE_N[di][ai],npass?"PASS":"FAIL",
                r.meanDeltaD,ARTICLE_DD[di][ai],dpassOne?"PASS":"DIFFERS");
        }
        System.out.printf("N cross-check: %d/%d pass%n",np,n);
        System.out.printf("DeltaD cross-check: %d/%d pass%n",dpass,dcount);
        System.out.println();
    }

    // ============================================================
    // MAIN
    // ============================================================
    public static void main(String[] args) {
        FAST_MODE = args.length > 0 && "--fast".equalsIgnoreCase(args[0]);
        System.out.println("EDBT 8-tuple: T+ = (V, E, Phi, Psi, Pi, gamma, Delta, Omega)");
        try { csv = new PrintWriter(new BufferedWriter(new FileWriter("edbt_pertrial.csv"))); 
              csv.println("depth,trial,alg,N,U,found,cost,deltaS_mean,deltaD,optCost"); } 
        catch (IOException e) { System.err.println("CSV open failed: " + e); }
        runExperiments();
        printArticleCrossCheck();
        if (csv != null) {
            for (String row : trialRows) csv.println(row);
            csv.close();
            System.out.println("Per-trial CSV written: edbt_pertrial.csv (" + trialRows.size() + " rows)");
        }
        if (!GraphicsEnvironment.isHeadless()) {
            SwingUtilities.invokeLater(() -> {
                JFrame f = new JFrame("EDBT Framework Results");
                f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                f.setSize(1400,700);
                UnifiedTreeExperimentL p = new UnifiedTreeExperimentL();
                p.setBackground(Color.WHITE);
                f.add(p);
                f.setLocationRelativeTo(null);
                f.setVisible(true);
            });
        }
    }
}
