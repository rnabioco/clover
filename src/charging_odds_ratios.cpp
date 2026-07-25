#include <cpp11.hpp>

#include <Rmath.h>

#include <cmath>
#include <vector>

using namespace cpp11;
using namespace cpp11::literals;

// Two-sided Fisher's exact test on a 2x2 table, as the sum of hypergeometric
// probabilities no larger than the observed one. Matches pairwise_fisher_exact.
static double fisher_two_sided(int a, int b, int c, int d) {
  int R1 = a + b;
  int R2 = c + d;
  int C1 = a + c;

  double p_obs = dhyper(static_cast<double>(a), static_cast<double>(R1),
                        static_cast<double>(R2), static_cast<double>(C1), 0);

  int a_min = (C1 - R2 > 0) ? C1 - R2 : 0;
  int a_max = (R1 < C1) ? R1 : C1;

  double p_value = 0.0;
  double rel_err = 1.0 + 1e-7;
  for (int t = a_min; t <= a_max; ++t) {
    double pt = dhyper(static_cast<double>(t), static_cast<double>(R1), static_cast<double>(R2),
                       static_cast<double>(C1), 0);
    if (pt <= p_obs * rel_err) {
      p_value += pt;
    }
  }

  return p_value > 1.0 ? 1.0 : p_value;
}

// Smallest two-sided p-value any outcome could produce for these margins.
//
// The extreme tables are the ones with the fewest ways of occurring, so their
// density bounds the attainable p-value from below. If even that bound exceeds
// the threshold the site cannot be called however the reads fall, and testing
// it only costs FDR budget.
static double min_attainable_p(int R1, int R2, int C1) {
  int a_min = (C1 - R2 > 0) ? C1 - R2 : 0;
  int a_max = (R1 < C1) ? R1 : C1;

  double p_lo = dhyper(static_cast<double>(a_min), static_cast<double>(R1),
                       static_cast<double>(R2), static_cast<double>(C1), 0);
  double p_hi = dhyper(static_cast<double>(a_max), static_cast<double>(R1),
                       static_cast<double>(R2), static_cast<double>(C1), 0);

  return p_lo < p_hi ? p_lo : p_hi;
}

// Accumulate modified-by-charged counts for every (reference, position) cell in
// a single pass over the long-form calls, then test each cell.
//
// The counts are the only state needed, so this never materializes a read x
// position matrix -- memory is O(n_refs * n_pos) regardless of read depth.
// Reads with no call at a position simply contribute no row there, so they are
// excluded from that position's table rather than counted as unmodified.
//
// Sites are pruned on their margins, never on the observed odds ratio: the
// margins are ancillary, so dropping sites by margin leaves the null
// distribution of the surviving p-values intact, while dropping them by effect
// size would enrich for small p and break the FDR correction downstream.
[[cpp11::register]]
writable::data_frame charging_odds_ratios_cpp(integers ref_idx, integers pos_idx,
                                              integers modified, integers charged,
                                              int n_refs, int n_pos, int min_reads,
                                              int min_margin, double max_p) {
  R_xlen_t n = ref_idx.size();

  std::size_t n_cells = static_cast<std::size_t>(n_refs) * static_cast<std::size_t>(n_pos);
  std::vector<int> n11(n_cells, 0), n10(n_cells, 0), n01(n_cells, 0), n00(n_cells, 0);

  for (R_xlen_t k = 0; k < n; ++k) {
    int r = ref_idx[k];
    int p = pos_idx[k];
    int m = modified[k];
    int ch = charged[k];

    if (r == NA_INTEGER || p == NA_INTEGER || m == NA_INTEGER || ch == NA_INTEGER) {
      continue;
    }

    // Inputs are 1-based factor codes.
    std::size_t cell =
        static_cast<std::size_t>(r - 1) * static_cast<std::size_t>(n_pos) + (p - 1);

    if (m) {
      if (ch) {
        n11[cell]++;
      } else {
        n10[cell]++;
      }
    } else {
      if (ch) {
        n01[cell]++;
      } else {
        n00[cell]++;
      }
    }
  }

  writable::integers ref_out;
  writable::integers pos_out;
  writable::doubles n11_out, n10_out, n01_out, n00_out, total_out;
  writable::doubles or_out, log_or_out, pval_out;

  for (std::size_t cell = 0; cell < n_cells; ++cell) {
    int a = n11[cell], b = n10[cell], c = n01[cell], d = n00[cell];
    int total = a + b + c + d;

    if (total < min_reads) {
      continue;
    }

    int R1 = a + b, R2 = c + d, C1 = a + c, C2 = b + d;

    // With an empty margin the odds ratio is undefined; a thin one carries no
    // power. min_margin >= 1 subsumes the empty case.
    if (R1 < min_margin || R2 < min_margin || C1 < min_margin || C2 < min_margin) {
      continue;
    }

    if (max_p < 1.0 && min_attainable_p(R1, R2, C1) > max_p) {
      continue;
    }

    double odds_ratio;
    if (a > 0 && b > 0 && c > 0 && d > 0) {
      odds_ratio = (static_cast<double>(a) * d) / (static_cast<double>(b) * c);
    } else {
      // Haldane correction for an empty cell.
      odds_ratio = (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5));
    }

    ref_out.push_back(static_cast<int>(cell / static_cast<std::size_t>(n_pos)) + 1);
    pos_out.push_back(static_cast<int>(cell % static_cast<std::size_t>(n_pos)) + 1);
    n11_out.push_back(a);
    n10_out.push_back(b);
    n01_out.push_back(c);
    n00_out.push_back(d);
    total_out.push_back(total);
    or_out.push_back(odds_ratio);
    log_or_out.push_back(std::log(odds_ratio));
    pval_out.push_back(fisher_two_sided(a, c, b, d));
  }

  return writable::data_frame({"ref_idx"_nm = ref_out, "pos_idx"_nm = pos_out,
                               "n11"_nm = n11_out, "n10"_nm = n10_out, "n01"_nm = n01_out,
                               "n00"_nm = n00_out, "total_obs"_nm = total_out,
                               "odds_ratio"_nm = or_out, "log_odds_ratio"_nm = log_or_out,
                               "p_value"_nm = pval_out});
}
