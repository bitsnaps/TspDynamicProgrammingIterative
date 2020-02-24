package project

/**
 * An implementation of the traveling salesman problem in Groovy using dynamic programming to improve
 * the time complexity from O(n!) to O(n^2 * 2^n).
 *
 * <p>Time Complexity: O(n^2 * 2^n) Space Complexity: O(n * 2^n)
 *
 * @author original Java code by William Fiset, william.alexandre.fiset@gmail.com
 */

class TspDynamicProgrammingIterative {

   final int N, start
   final double[][] distance
   def tour = []
   double minTourCost = Double.POSITIVE_INFINITY
   def ranSolver = false

  TspDynamicProgrammingIterative(double[][] distance) {
    this(0, distance)
  }

  TspDynamicProgrammingIterative(int start, double[][] distance) {
    N = distance.length

      def throwException = { def m ->
          throw new IllegalStateException(m)
          
      }
    if (N <= 2) throwException("N <= 2 not yet supported.")
    if (N != distance[0].length) 
      throwException("Matrix must be square (n x n)")
    if (start < 0 || start >= N) 
      throwException("Invalid start node.")
    if (N > 32)
      throwException(
          "Matrix too large! A matrix that size for the DP TSP problem with a time complexity of"
              + "O(n^2*2^n) requires way too much computation for any modern home computer to handle")

    this.start = start
    this.distance = distance
  }

  // Returns the optimal tour for the traveling salesman problem.
  List<Integer> getTour() {
    if (!ranSolver) solve()
    tour
  }

  // Returns the minimal tour cost.
  double getTourCost() {
    if (!ranSolver) solve()
    minTourCost
  }

  // Solves the traveling salesman problem and caches solution.
  void solve() {

    if (ranSolver) return

    final int END_STATE = (1 << N) - 1
    Double[][] memo = new Double[N][1 << N]

    // Add all outgoing edges from the starting node to memo table.
    for (int end = 0; end < N; end++) {
      if (end == start) continue
      memo[end][(1 << start) | (1 << end)] = distance[start][end]
    }

    for (int r = 3; r <= N; r++) {
      for (int subset : combinations(r, N)) {
        if (notIn(start, subset)) continue
        for (int next = 0; next < N; next++) {
          if (next == start || notIn(next, subset)) continue
          int subsetWithoutNext = subset ^ (1 << next)
          double minDist = Double.POSITIVE_INFINITY
          for (int end = 0; end < N; end++) {
            if (end == start || end == next || notIn(end, subset)) continue
            double newDistance = memo[end][subsetWithoutNext] + distance[end][next];
            if (newDistance < minDist) {
              minDist = newDistance;
            }
          }
          memo[next][subset] = minDist
        }
      }
    }

    // Connect tour back to starting node and minimize cost.
    for (int i = 0; i < N; i++) {
      if (i == start) continue;
      double tourCost = memo[i][END_STATE] + distance[i][start]
      if (tourCost < minTourCost) {
        minTourCost = tourCost
      }
    }

    int lastIndex = start
    int state = END_STATE
    tour << start

    // Reconstruct TSP path from memo table.
    for (int i = 1; i < N; i++) {

      int index = -1
      for (int j = 0; j < N; j++) {
        if (j == start || notIn(j, state)) continue
        if (index == -1) index = j
        double prevDist = memo[index][state] + distance[index][lastIndex]
        double newDist = memo[j][state] + distance[j][lastIndex]
        if (newDist < prevDist) {
          index = j
        }
      }

      tour << index
      state = state ^ (1 << index)
      lastIndex = index
    }

    tour << start
    tour.reverse()

    ranSolver = true
  }

   static boolean notIn(int elem, int subset) {
    ((1 << elem) & subset) == 0
  }

  // This method generates all bit sets of size n where r bits
  // are set to one. The result is returned as a list of integer masks.
  public static List<Integer> combinations(int r, int n) {
    List<Integer> subsets = []
    combinations(0, 0, r, n, subsets)
    subsets
  }

  // To find all the combinations of size r we need to recurse until we have
  // selected r elements (aka r = 0), otherwise if r != 0 then we still need to select
  // an element which is found after the position of our last selected element
  static void combinations(int set, int at, int r, int n, List<Integer> subsets) {

    // Return early if there are more elements left to select than what is available.
    int elementsLeftToPick = n - at;
    if (elementsLeftToPick < r) return;

    // We selected 'r' elements so we found a valid subset!
    if (r == 0) {
      subsets << set
    } else {
      for (int i = at; i < n; i++) {
        // Try including this element
        set ^= (1 << i)

        combinations(set, i + 1, r - 1, n, subsets)

        // Backtrack and try the instance where we did not include this element
        set ^= (1 << i)
      }
    }
  }

  static void main(String[] args) {
    // Create adjacency matrix
    int n = 6
    double[][] distanceMatrix = new double[n][n]
    //for (double[] row : distanceMatrix)
      distanceMatrix.each { def row -> Arrays.fill(row, 10000) }
    distanceMatrix[5][0] = 10
    distanceMatrix[1][5] = 12
    distanceMatrix[4][1] = 2
    distanceMatrix[2][4] = 4
    distanceMatrix[3][2] = 6
    distanceMatrix[0][3] = 8

    int startNode = 0
    TspDynamicProgrammingIterative solver = new TspDynamicProgrammingIterative(startNode, distanceMatrix)

    // Prints: [0, 3, 2, 4, 1, 5, 0]
    def tour = solver.tour
    println("Tour: " + tour)
    assert tour == [0, 5, 1, 4, 2, 3, 0]
    

    // Print: 42.0
    def tourCost = solver.tourCost
    println("Tour cost: " + tourCost)
    assert tourCost == 42.0
    
  }

    
}
