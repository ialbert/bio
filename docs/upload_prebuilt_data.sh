set -uex

BUCKET_NAME="gs://biostore-bucket-001"

# This directory will sync the prebuilt data
PRE="~/.bio/pre"

mkdir -p $PRE

# The files that need to be synchronized to the bucket.
ln -sf ~/.bio/taxdb.json $PRE

ln -sf ~/.bio/taxdb.sqlite $PRE

time gsutil rsync $PRE $BUCKET_NAME
