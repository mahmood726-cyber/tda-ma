## REVIEW CLEAN
## Code Review: tda_ma.html
### Date: 2026-04-03
### Summary: 1 P0 FIXED, 1 P1, 3 P2

#### P0 -- Critical (FIXED)
- **[FIXED] P0-1** [Security]: CSV export (`exportCSV`) writes study names and column values directly into CSV cells without formula injection protection. Custom CSV input allows arbitrary study names starting with `=`, `+`, `@`, `\t`, `\r` -- these can trigger formula execution when opened in Excel/Sheets. Fix: add `sanitizeCSVCell()` to prepend `'` to cells starting with dangerous characters (per rule: NOT `-`).

#### P1 -- Important
- **P1-1** [Stats]: `standardize()` uses population SD (`ssq / k`) instead of sample SD (`ssq / (k-1)`). For z-score normalization of study data, population SD is defensible (consistent with the intended normalization semantics), but differs from R's `scale()` which uses sample SD. The difference is negligible for k >= 10 (BCG k=13, Magnesium k=16). Left as-is since it's consistent with the documented behavior ("Population SD (consistent with z-score normalization)") and the PCA covariance correctly uses `k-1`.

#### P2 -- Minor
- **P2-1** [Stats]: DerSimonian-Laird pooling (`dlPooled`) correctly implements inverse-variance weighting with `tau2 = max(0, (Q - (k-1)) / C)`. The I-squared formula `max(0, (Q - (k-1)) / Q) * 100` is correct.
- **P2-2** [Stats]: Persistent homology H1 via boundary matrix reduction is correctly implemented with Z/2 coefficients. The square test case validates: H1 born at eps=1, dies at sqrt(2).
- **P2-3** [Code Quality]: `escapeHtml` correctly escapes all 5 characters, verified in test 12. Tooltip innerHTML uses `escapeHtml()` for study names (line 2830). Dimension checkbox creation uses `escapeHtml()` for column names (line 2159).

#### Structural Checks
- Div balance: 45 open, 45 close -- balanced
- `</script>` literal: only at line 3617 (actual closing tag) -- safe
- `</html>` closing tag: present at line 3619
- PRNG: uses `xoshiro128**` with seed -- deterministic
- Blob URLs: `URL.revokeObjectURL()` called after download -- no leak
- Modal close: keydown listener properly added/removed in `initTestModal()`

#### Test Results: 25/25 pass (in-browser suite)
