#!/bin/bash
# Mirror human CELLxGENE Census 2025-11-08 SOMA store (~1.2 TB). Re-run to resume.
aws configure set default.s3.max_concurrent_requests 32
aws configure set default.s3.multipart_chunksize 64MB
aws s3 sync --no-sign-request --region us-west-2 --exclude "census_data/mus_musculus/*" s3://cellxgene-census-public-us-west-2/cell-census/2025-11-08/soma/ /mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma/
echo SYNC_EXIT=$?
