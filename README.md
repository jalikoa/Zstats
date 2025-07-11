# ZStats – Comprehensive PHP Statistics Library

## Features

- Z-score calculator
- Confidence interval builder
- T-distribution & Chi-square approximations
- HTML/CSV/JSON exports
- Unit tested

## Install via Composer

```bash
composer require jalsoft/zstats
```
Usage
```php
use ZStats\ZStats;

$ci = ZStats::confidenceInterval(3.3, 0.2, 64, 0.98);
print_r($ci);
```
