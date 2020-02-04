#!/bin/bash
conda activate HeuDiConv

docker run --rm -it -v /Volumes/DATA:/base nipy/heudiconv:latest -d /base/sourcedata/sub-{subject}/ses-{session}/*/*.dcm -o /base/out -f convertall -s 48 -ss 01 -c none --overwrite
docker run --rm -it -v /Volumes/DATA:/base nipy/heudiconv:latest -d /base/sourcedata/sub-{subject}/ses-{session}/*/*.dcm -o /base/out -f /base/code/heuristic.py -s 48 -ss 01 -c dcm2niix -b --overwrite

docker run --rm -it -v /Volumes/DATA:/base nipy/heudiconv:latest -d /base/sourcedata/sub-{subject}/ses-{session}/*/*.dcm -o /base/out -f convertall -s 50 -ss 01 -c none --overwrite
docker run --rm -it -v /Volumes/DATA:/base nipy/heudiconv:latest -d /base/sourcedata/sub-{subject}/ses-{session}/*/*.dcm -o /base/out -f /base/code/heuristic.py -s 50 -ss 01 -c dcm2niix -b --overwrite

docker run --rm -it -v /Volumes/DATA:/base nipy/heudiconv:latest -d /base/sourcedata/sub-{subject}/ses-{session}/*/*.dcm -o /base/out -f convertall -s 82 -ss 01 -c none --overwrite
docker run --rm -it -v /Volumes/DATA:/base nipy/heudiconv:latest -d /base/sourcedata/sub-{subject}/ses-{session}/*/*.dcm -o /base/out -f /base/code/heuristic.py -s 82 -ss 01 -c dcm2niix -b --overwrite

docker run --rm -it -v /Volumes/DATA:/base nipy/heudiconv:latest -d /base/sourcedata/sub-{subject}/ses-{session}/*/*.dcm -o /base/out -f convertall -s 96 -ss 01 -c none --overwrite
docker run --rm -it -v /Volumes/DATA:/base nipy/heudiconv:latest -d /base/sourcedata/sub-{subject}/ses-{session}/*/*.dcm -o /base/out -f /base/code/heuristic.py -s 96 -ss 01 -c dcm2niix -b --overwrite
