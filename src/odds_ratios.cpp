#include <cpp11.hpp>

#include <Rmath.h>

using namespace cpp11;
using namespace cpp11::literals;

[[cpp11::register]]
writable::data_frame pairwise_fisher_exact(integers_matrix<> mat) {
  int nr = mat.nrow();
  int nc = mat.ncol();

  writable::integers pos1_out;
  writable::integers pos2_out;
  writable::doubles or_out;
  writable::doubles log_or_out;
  writable::doubles pval_out;
  writable::integers total_out;

  for (int i = 0; i < nc - 1; ++i) {
    for (int j = i + 1; j < nc; ++j) {
      int a = 0, b = 0, c = 0, d = 0;
      for (int k = 0; k < nr; ++k) {
        int vi = mat(k, i);
        int vj = mat(k, j);
        if (vi == 1 && vj == 1) {
          a++;
        } else if (vi == 1 && vj == 0) {
          b++;
        } else if (vi == 0 && vj == 1) {
          c++;
        } else {
          d++;
        }
      }

      // Skip if any marginal is zero (no variation in one direction)
      int R1 = a + b;
      int R2 = c + d;
      int C1 = a + c;
      int C2 = b + d;

      if (R1 == 0 || R2 == 0 || C1 == 0 || C2 == 0) {
        continue;
      }

      // Sample odds ratio with Haldane correction for zero cells
      double odds_ratio;
      if (a > 0 && b > 0 && c > 0 && d > 0) {
        odds_ratio = (static_cast<double>(a) * d) / (static_cast<double>(b) * c);
      } else {
        odds_ratio = (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5));
      }

      // Fisher's exact test (two-sided, sum of small probabilities)
      // a ~ Hypergeometric(R1, R2, C1)
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

      if (p_value > 1.0)
        p_value = 1.0;

      pos1_out.push_back(i + 1);
      pos2_out.push_back(j + 1);
      or_out.push_back(odds_ratio);
      log_or_out.push_back(std::log(odds_ratio));
      pval_out.push_back(p_value);
      total_out.push_back(nr);
    }
  }

  return writable::data_frame({"pos1"_nm = pos1_out, "pos2"_nm = pos2_out, "odds_ratio"_nm = or_out,
                               "log_odds_ratio"_nm = log_or_out, "p_value"_nm = pval_out,
                               "total_obs"_nm = total_out});
}
