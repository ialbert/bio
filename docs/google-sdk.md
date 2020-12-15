
From:

* https://cloud.google.com/sdk/docs/downloads-interactive


        # Install the Google Cloud.
        curl https://sdk.cloud.google.com > install.sh
        
        bash install.sh --disable-prompts
    
        # Install GS Util (google cloud file browser)
        pip install gsutil
        
        # Initialize the Google Cloud
        gcloud init
       
        # Bucket 
        BUCKET_NAME=gs://biostore-bucket-001
        FILE_NAME=foo.txt
        
        # Make the storage bucket Public
        gsutil iam ch allUsers:objectViewer $BUCKET_NAME
        
        # Upload
        gsutil cp $FILE_NAME $BUCKET_NAME
        
        # Download    
        wget "https://storage.googleapis.com/storage/v1/b/biostore-bucket-001/o/$FILE_NAME"
        
        # Delete 
        gsutil rm $BUCKET_NAME/$FILE_NAME
      
  
Disable annoying upload message:

Change this section in `~/.boto`:
      
    [GSUtil]
    parallel_composite_upload_threshold = 550M

 