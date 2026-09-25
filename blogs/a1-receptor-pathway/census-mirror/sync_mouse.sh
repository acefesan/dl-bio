#!/bin/bash
# Mirror mouse CELLxGENE Census 2025-11-08 SOMA store (~281 GB). Re-run to resume.
aws s3 sync --no-sign-request --region us-west-2 s3://cellxgene-census-public-us-west-2/cell-census/2025-11-08/soma/census_data/mus_musculus/ /mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma/census_data/mus_musculus/
echo SYNC_EXIT=$?
