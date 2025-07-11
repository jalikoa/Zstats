<?php
require_once '../src/ZStats.php';

use PHPUnit\Framework\TestCase;
use ZStats\ZStats;

class ZStatsTest extends TestCase
{
    public function testZScore()
    {
        $this->assertEqualsWithDelta(1.0, ZStats::zScore(110, 100, 10), 0.001);
    }

    public function testMarginOfError()
    {
        $this->assertEqualsWithDelta(1.96, ZStats::marginOfError(10, 100, 1.96), 0.001);
    }

    public function testConfidenceInterval()
    {
        $ci = ZStats::confidenceInterval(50, 10, 100, 0.95);
        $this->assertGreaterThan(48, $ci['lower']);
        $this->assertLessThan(52, $ci['upper']);
    }

    public function testTcdf()
    {
        $this->assertGreaterThan(0.5, ZStats::tDistributionCDF(1, 10));
    }

    public function testChiSquareCDF()
    {
        $this->assertGreaterThan(0.5, ZStats::chiSquareCDF(5, 3));
    }

    public function testZTableExport()
    {
        $table = ZStats::generateZTable();
        $this->assertArrayHasKey('0', $table);$found = false;
        foreach ($table as $z => $prob) {
            if (abs($z - 1.0) < 0.01) { // Look for any key within 0.01 of 1.0
                $this->assertGreaterThan(0.5, $prob);
                $found = true;
                break;
            }
        }
        $this->assertTrue($found, "No Z-value close to 1.0 found in Z-table");
    }
}