<?php

require_once __DIR__ . '/../vendor/autoload.php';

use ZStats\ZStats;

echo "📊 Example Usage of ZStats Library\n";
echo "==================================\n\n";

// 🔹 Basic Z-score calculation
$x = 85;
$mean = 80;
$stdDev = 5;
$zScore = ZStats::zScore($x, $mean, $stdDev);
echo "🔹 Z-Score for x=$x, mean=$mean, stdDev=$stdDev => $zScore\n\n";

// 🔹 Confidence Interval
$sampleMean = 78;
$populationStdDev = 6;
$sampleSize = 100;
$confidenceLevel = 0.95;

$ci = ZStats::confidenceInterval($sampleMean, $populationStdDev, $sampleSize, $confidenceLevel);
echo "🔹 Confidence Interval ($confidenceLevel):\n";
echo "   Lower: " . $ci['lower'] . "\n";
echo "   Upper: " . $ci['upper'] . "\n";
echo "   Margin of Error: " . $ci['margin_of_error'] . "\n\n";

// 🔹 Critical Z-value for 95% confidence
$alpha = 1 - $confidenceLevel;
$zCritical = ZStats::criticalZValue($alpha);
echo "🔹 Critical Z-value for alpha=$alpha => $zCritical\n\n";

// 🔹 P-value from Z-score
$pValue = ZStats::pValueFromZScore($zScore);
echo "🔹 P-value for Z-score $zScore => $pValue\n\n";

// 🔹 T-distribution CDF
$t = 2.0;
$df = 10;
$tCdf = ZStats::tDistributionCDF($t, $df);
echo "🔹 T-distribution CDF for t=$t, df=$df => $tCdf\n\n";

// 🔹 Chi-square distribution CDF
$chi2 = 5.0;
$chiDf = 3;
$chiCdf = ZStats::chiSquareCDF($chi2, $chiDf);
echo "🔹 Chi-square CDF for x²=$chi2, df=$chiDf => $chiCdf\n\n";