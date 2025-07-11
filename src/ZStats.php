<?php

namespace ZStats;

/**
 * A comprehensive statistical library with support for normal, t, chi-square distributions,
 * confidence intervals, hypothesis testing, and data exporting.
 */
class ZStats
{
    // --- Normal Distribution Functions ---

    /**
     * Compute the error function erf(x)
     */
    public static function erf(float $x): float
    {
        $sign = ($x < 0) ? -1 : 1;
        $x = abs($x);
        $t = 1 / (1 + 0.5 * $x);
        $tau = $t * exp(-$x * $x - 1.26551223 +
            $t * (1.00002368 +
            $t * (0.37409196 +
            $t * (0.09678418 +
            $t * (-0.18628806 +
            $t * (0.27886807 +
            $t * (-1.13520398 +
            $t * (1.48851587 +
            $t * (-0.82215223 +
            $t * 0.17087277)))))))));

        return $sign * (1 - $tau);
    }

    /**
     * Standard normal cumulative distribution function
     */
    public static function normalCDF(float $z): float
    {
        return 0.5 * (1 + self::erf($z / sqrt(2)));
    }

    /**
     * Inverse standard normal CDF (probit function)
     */
    public static function inverseNormalCDF(float $p, float $mu = 0, float $sigma = 1): float
    {
        if ($p <= 0 || $p >= 1) {
            throw new \InvalidArgumentException("Probability must be strictly between 0 and 1.");
        }

        $a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
              1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];

        $b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
              6.680131188771972e+01, -1.328068155288572e+01];

        $c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
              -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];

        $d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
              3.754408661907416e+00];

        $plow = 0.02425;
        $phigh = 1 - $plow;

        if ($p < $plow) {
            $q = sqrt(-2 * log($p));
            $z = ((((($c[0] * $q + $c[1]) * $q + $c[2]) * $q + $c[3]) * $q + $c[4]) * $q + $c[5]) /
                  ((((($d[0] * $q + $d[1]) * $q + $d[2]) * $q + $d[3]) * $q + 1));
        } elseif ($p > $phigh) {
            $q = sqrt(-2 * log(1 - $p));
            $z = -(((((($c[0] * $q + $c[1]) * $q + $c[2]) * $q + $c[3]) * $q + $c[4]) * $q + $c[5]) /
                   ((((($d[0] * $q + $d[1]) * $q + $d[2]) * $q + $d[3]) * $q + 1)));
        } else {
            $q = $p - 0.5;
            $r = $q * $q;
            $z = ((((($a[0] * $r + $a[1]) * $r + $a[2]) * $r + $a[3]) * $r + $a[4]) * $r + $a[5]) * $q /
                 ((((($b[0] * $r + $b[1]) * $r + $b[2]) * $r + $b[3]) * $r + $b[4]) * $r + 1);
        }

        return $mu + $sigma * $z;
    }

    // --- New Additions Below ---

    /**
     * Basic z-score calculation
     */
    public static function zScore(float $x, float $mean, float $stdDev): float
    {
        if ($stdDev <= 0) {
            throw new \InvalidArgumentException("Standard deviation must be positive.");
        }
        return ($x - $mean) / $stdDev;
    }

    /**
     * Calculate margin of error
     */
    public static function marginOfError(float $stdDev, int $n, float $z): float
    {
        if ($n <= 0) {
            throw new \InvalidArgumentException("Sample size must be positive.");
        }
        return $z * ($stdDev / sqrt($n));
    }

    /**
     * Confidence interval for population mean
     */
    public static function confidenceInterval(float $mean, float $stdDev, int $n, float $confidenceLevel): array
    {
        if ($confidenceLevel <= 0 || $confidenceLevel >= 1) {
            throw new \InvalidArgumentException("Confidence level must be between 0 and 1.");
        }

        $alpha = 1 - $confidenceLevel;
        $z = self::inverseNormalCDF(1 - $alpha / 2);

        $error = self::marginOfError($stdDev, $n, $z);
        return [
            'lower' => $mean - $error,
            'upper' => $mean + $error,
            'margin_of_error' => $error,
            'z_score' => $z
        ];
    }

    /**
     * Approximate t-distribution CDF using series expansion
     */
    public static function tDistributionCDF(float $t, int $df): float
    {
        if ($df <= 0) {
            throw new \InvalidArgumentException("Degrees of freedom must be positive.");
        }

        // Use regularized incomplete beta function
        $x = $t * $t / ($df + $t * $t);
        $a = $df / 2;
        $b = 0.5;

        $beta_func = self::incompleteBetaFunction($x, $a, $b);

        return $t < 0 ? 0.5 * $beta_func : 1 - 0.5 * $beta_func;
    }

    private static function lnBeta(float $a, float $b): float
    {
        return self::lgamma($a) + self::lgamma($b) - self::lgamma($a + $b);
    }

    private static function incompleteBetaFunction(float $x, float $a, float $b): float
    {
        if ($x < 0 || $x > 1) return 0;
        if ($a <= 0 || $b <= 0) return 0;

        $BTOL = 1e-11;
        $MAXIT = 1000;

        if ($x == 0 || $x == 1) {
            return $x;
        }

        $beta = exp(self::lnBeta($a, $b));
        $bt = $beta * pow($x, $a) * pow(1 - $x, $b);

        $p = 1;
        $q = 1;
        $m = 0;
        $m2 = 0;

        $aa = 0;
        $c = 1;
        $d = 1 - (1 + $a + $b) * $x / ($a + 1);

        if (abs($d) < 1e-30) $d = 1e-30;
        $d = 1 / $d;
        $h = $d;

        for ($i = 1; $i <= $MAXIT; $i++) {
            $m = $i;
            $m2 = 2 * ($i - 1);
            $aa = $m2 + 1;

            // Even step
            $d = 1 + ($m + $b) * ($aa - $a) * $x / (($aa) * ($aa + 1) / ($m + $a));
            if (abs($d) < 1e-30) $d = 1e-30;
            $d = 1 / $d;
            $h *= $d;
            $c = 1 + $m * ($aa - $a) * $x / (($aa + 1) * ($aa) / ($m + $a));
            if (abs($c) < 1e-30) $c = 1e-30;
            $h /= $c;
            if (abs($h - 1) < $BTOL) break;
        }

        return $bt * $h / $a;
    }
    /**
     * Natural logarithm of the Gamma function (lgamma)
     * Approximation based on the Lanczos algorithm.
     *
     * @param float $x Input value > 0
     * @return float ln(Gamma(x))
     */
    public static function lgamma(float $x): float
    {
        if ($x <= 0.0) {
            throw new \InvalidArgumentException("Argument must be positive.");
        }

        static $g = 7;
        static $c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                    771.32342877765313, -176.61502916214059, 12.507341404690872,
                    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];

        $xx = $x;
        $t = $xx + $g + 0.5;

        $sum = $c[0];
        for ($i = 1; $i < count($c); $i++) {
            $sum += $c[$i] / ($xx + $i);
        }

        return (0.5 * M_LNPI) + \log($sum) - $xx - ($xx - 0.5) * \log($t);
    }
    /**
     * Chi-square CDF approximation using gamma function
     */
    public static function chiSquareCDF(float $x, int $df): float
    {
        if ($x < 0 || $df <= 0) {
            throw new \InvalidArgumentException("Invalid input for chi-square distribution.");
        }

        $gamma = function ($z) {
            return exp(-$z) * pow($z, $z - 0.5) * sqrt(2 * M_PI);
        };

        $g = 0;
        $sum = 0;
        for ($k = 0; $k <= $df / 2 - 1; $k++) {
            $sum += pow($x / 2, $k) / self::factorial($k);
        }
        $g = 1 - exp(-$x / 2) * $sum;

        return $g;
    }

    private static function factorial(int $n): float
    {
        $result = 1;
        for ($i = 2; $i <= $n; $i++) {
            $result *= $i;
        }
        return $result;
    }

    /**
     * Format Z-table as HTML
     */
    public static function formatZTableAsHtml(array $zTable): string
    {
        $html = "<table border='1'><thead><tr><th>Z</th><th>Cumulative Probability</th></tr></thead><tbody>";
        foreach ($zTable as $z => $prob) {
            $html .= "<tr><td>$z</td><td>" . round($prob, 5) . "</td></tr>";
        }
        $html .= "</tbody></table>";
        return $html;
    }

    /**
     * Export Z-table to JSON
     */
    public static function exportZTableToJson(array $zTable, string $filePath): void
    {
        file_put_contents($filePath, json_encode($zTable, JSON_PRETTY_PRINT));
    }

    /**
     * Export Z-table to CSV
     */
    public static function exportZTableToCsv(array $zTable, string $filePath): void
    {
        $fp = fopen($filePath, 'w');
        fputcsv($fp, ['Z', 'Cumulative Probability']);
        foreach ($zTable as $z => $prob) {
            fputcsv($fp, [$z, $prob]);
        }
        fclose($fp);
    }

    /**
     * Generate Z-table
     */
    public static function generateZTable(float $min = -3.9, float $max = 3.9, float $step = 0.01): array
    {
        $table = [];
        for ($z = $min; $z <= $max; $z += $step) {
            $roundedZ = round($z, 2);
            $table[$roundedZ] = round(self::normalCDF($roundedZ), 5);
        }
        return $table;
    }
}