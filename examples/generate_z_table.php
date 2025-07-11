<?php

// Make sure autoloader is working
require_once __DIR__ . '/../vendor/autoload.php';

use ZStats\ZStats;

// Generate Z-table from -3.9 to 3.9 with step size of 0.01
$zTable = ZStats::generateZTable(-3.9, 3.9, 0.01);

// Export to JSON
ZStats::exportZTableToJson($zTable, __DIR__ . '/../z_table.json');
echo "✅ Z-table exported to z_table.json\n";

// Export to CSV
ZStats::exportZTableToCsv($zTable, __DIR__ . '/../z_table.csv');
echo "✅ Z-table exported to z_table.csv\n";