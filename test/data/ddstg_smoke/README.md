# DDSTG Smoke Fixture

This directory contains the repo-local synthetic fixture used by
`scripts/smoke-run.sh`.

Files:

- short default smoke fixture:
  - `j1141-6545_ddstg_gr_pdfb1_white_noise.par`
  - `j1141-6545_ddstg_gr_pdfb1_white_noise.tim`
- larger multi-backend synthetic fixture:
  - `j1141-6545_ddstg_gr_white_noise.par`
  - `j1141-6545_ddstg_gr_white_noise.tim`

Fixture intent:

- validate DDSTG `Tempo2` wiring in the rebuilt development container
- provide a repository-local smoke dataset that does not depend on private
  host-mounted research directories
- keep both a short default fixture and a broader synthetic fixture available
  for future tests

Provenance:

- derived from synthetic white-noise data based on a real system
- approved for repository use as test data
- not raw observational data

Sanitization:

- the larger multi-backend `.tim` fixture had its first-column filesystem paths
  reduced to basenames so internal machine paths are not retained in the
  repository
- the short `PDFB1` `.tim` fixture already used basename-style file identifiers

Notes:

- `scripts/smoke-run.sh` defaults to the shorter `PDFB1` fixture
- the larger multi-backend fixture remains available for broader DDSTG testing
- smoke-run output is written under `DDSTG_SMOKE/`, which is ignored by git
